"""Offline SCC credential-prompt and authenticated-download regression tests."""
import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch
import warnings

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from atlas_actris.utils import get_scc_config as scc
from atlas_actris import __get_config_file_from_scc_hoi__ as cli


class SccAuthenticationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.url = scc.SCC_EXPORT_BASE_URL + "/713"

    @staticmethod
    def response(status=200, chunks=(b"configuration data\n",), content_type="text/csv"):
        response = MagicMock()
        response.status_code = status
        response.headers = {"Content-Type": content_type}
        response.iter_content.return_value = iter(chunks)
        response.__enter__.return_value = response
        return response

    def test_hidden_prompt_keeps_password_whitespace(self):
        with patch("builtins.input", return_value=" test-user "):
            with patch.object(scc.getpass, "getpass", return_value=" test-password ") as hidden:
                self.assertEqual(scc._prompt_scc_credentials(), ("test-user", " test-password "))
        hidden.assert_called_once()
        self.assertEqual(list(self.folder.iterdir()), [])

    def test_visible_fallback_is_not_permitted(self):
        fallback = MagicMock()
        def unsafe_input(*args, **kwargs):
            warnings.warn("Cannot disable echo", scc.getpass.GetPassWarning)
            return fallback()
        with patch("builtins.input", return_value="test-user"):
            with patch.object(scc.getpass, "getpass", side_effect=unsafe_input):
                with self.assertRaisesRegex(scc.SCCAuthenticationError, "Hidden password input"):
                    scc._prompt_scc_credentials()
        fallback.assert_not_called()

    def test_unavailable_and_empty_credentials(self):
        with patch("builtins.input", side_effect=EOFError):
            with self.assertRaises(scc.SCCAuthenticationError):
                scc._prompt_scc_credentials()
        for username, password in (("", "test-password"), ("bad:name", "test-password"), ("test-user", "")):
            with self.subTest(username=username), patch("builtins.input", return_value=username):
                with patch.object(scc.getpass, "getpass", return_value=password):
                    with self.assertRaises(scc.SCCAuthenticationError):
                        scc._prompt_scc_credentials()

    def test_each_download_prompts_without_reusing_credentials(self):
        with patch.object(scc, "_prompt_scc_credentials", side_effect=[("user-one", "pass-one"), ("user-two", "pass-two")]) as prompt:
            with patch.object(scc, "download_with_auth", return_value="test.csv") as download:
                for _ in range(2):
                    scc._get_scc_file(str(self.folder), "713")
        self.assertEqual(prompt.call_count, 2)
        self.assertEqual([c.kwargs["username"] for c in download.call_args_list], ["user-one", "user-two"])
        self.assertEqual([c.kwargs["password"] for c in download.call_args_list], ["pass-one", "pass-two"])

    def test_no_prompt_for_disabled_export_or_existing_local_configuration(self):
        config = self.folder / "config.ini"
        config.touch()
        with patch.object(scc, "_prompt_scc_credentials") as prompt:
            with patch.object(scc.requests, "get") as network:
                for mode in ("0", "2"):
                    scc.export_scc_config("713", str(config), mode, str(self.folder))
        prompt.assert_not_called()
        network.assert_not_called()

    def test_configuration_id_checked_before_prompt(self):
        with patch.object(scc, "_prompt_scc_credentials") as prompt:
            for value in ("../713", "713?query=x", "-1", "0", "abc"):
                with self.assertRaises(ValueError):
                    scc._get_scc_file(str(self.folder), value)
        prompt.assert_not_called()
        self.assertEqual(scc._normalise_configuration_id("000713"), "713")

    def test_download_uses_https_timeout_and_atomic_write(self):
        response = self.response()
        with patch.object(scc.requests, "get", return_value=response) as request:
            path = scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password", "result.csv")
        self.assertEqual(Path(path).read_bytes(), b"configuration data\n")
        self.assertEqual(list(self.folder.iterdir()), [Path(path)])
        kwargs = request.call_args.kwargs
        self.assertEqual(kwargs["auth"], ("test-user", "test-password"))
        self.assertEqual(kwargs["timeout"], (10, 60))
        self.assertIs(kwargs["verify"], True)
        self.assertIs(kwargs["allow_redirects"], False)

    def test_unexpected_urls_refused_before_sending_credentials(self):
        with patch.object(scc.requests, "get") as request:
            for url in (self.url.replace("https:", "http:"),
                        self.url.replace("scc.imaa.cnr.it", "example.invalid"),
                        self.url + "?redirect=other", self.url.replace("713", "../other")):
                with self.assertRaises(scc.SCCDownloadError):
                    scc.download_with_auth(url, str(self.folder), "test-user", "test-password")
        request.assert_not_called()

    def test_auth_rejection_and_redirect_do_not_write_files(self):
        for status in (401, 403, 302, 500):
            with self.subTest(status=status):
                with patch.object(scc.requests, "get", return_value=self.response(status)) as request:
                    with self.assertRaises(scc.SCCDownloadError):
                        scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password", "result.csv")
                request.assert_called_once()
                self.assertEqual(list(self.folder.iterdir()), [])

    def test_canonical_redirects_keep_credentials_on_same_https_export(self):
        target = self.url + "/"
        for status in (301, 302, 303, 307, 308):
            with self.subTest(status=status):
                redirect = self.response(status)
                redirect.headers["Location"] = target.replace("https://scc.imaa.cnr.it", "")
                with patch.object(scc.requests, "get", side_effect=[redirect, self.response()]) as request:
                    result = scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password",
                                                   "{}.csv".format(status))
                self.assertEqual(Path(result).read_bytes(), b"configuration data\n")
                self.assertEqual([call.args[0] for call in request.call_args_list], [self.url, target])
                for call in request.call_args_list:
                    self.assertEqual(call.kwargs["auth"], ("test-user", "test-password"))
                    self.assertIs(call.kwargs["verify"], True)
                    self.assertIs(call.kwargs["allow_redirects"], False)
                    self.assertEqual(call.kwargs["timeout"], (10, 60))

    def test_explicit_default_https_port_is_same_origin(self):
        target = self.url.replace("scc.imaa.cnr.it", "scc.imaa.cnr.it:443") + "/"
        redirect = self.response(301)
        redirect.headers["Location"] = target
        with patch.object(scc.requests, "get", side_effect=[redirect, self.response()]) as request:
            scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password", "result.csv")
        self.assertEqual(request.call_count, 2)
        self.assertEqual(request.call_args.args[0], target)

    def test_unsafe_redirect_targets_never_receive_credentials(self):
        targets = (
            self.url.replace("https:", "http:"),
            self.url.replace("scc.imaa.cnr.it", "example.invalid"),
            "//example.invalid/export/713",
            self.url.replace("713", "714"),
            self.url.replace("scc.imaa.cnr.it", "scc.imaa.cnr.it:444"),
            self.url.replace("scc.imaa.cnr.it", "scc.imaa.cnr.it:0"),
            self.url.replace("https://", "https://test-user:test-password@"),
            self.url + "?token=test-password",
            self.url + "#test-password",
            "https://scc.imaa.cnr.it/admin/login/?next=/private&token=test-password",
            self.url + "/\n",
        )
        for target in targets:
            with self.subTest(target=target):
                redirect = self.response(302)
                redirect.headers["Location"] = target
                with patch.object(scc.requests, "get", return_value=redirect) as request:
                    with self.assertRaises(scc.SCCDownloadError) as caught:
                        scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password", "result.csv")
                request.assert_called_once()
                self.assertNotIn("test-password", str(caught.exception))
                self.assertEqual(list(self.folder.iterdir()), [])

    def test_unsafe_later_redirect_preserves_previous_file(self):
        final = self.folder / "result.csv"
        final.write_bytes(b"old valid data")
        first, second = self.response(301), self.response(302)
        first.headers["Location"] = self.url + "/"
        second.headers["Location"] = "https://example.invalid/private"
        with patch.object(scc.requests, "get", side_effect=[first, second]) as request:
            with self.assertRaises(scc.SCCDownloadError):
                scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password", final.name)
        self.assertEqual(request.call_count, 2)
        self.assertEqual(final.read_bytes(), b"old valid data")
        self.assertEqual(list(self.folder.iterdir()), [final])

    def test_redirect_loop_has_a_fixed_limit(self):
        redirect = self.response(302)
        redirect.headers["Location"] = self.url
        with patch.object(scc.requests, "get", return_value=redirect) as request:
            with self.assertRaisesRegex(scc.SCCDownloadError, "redirect limit"):
                scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password", "result.csv")
        self.assertEqual(request.call_count, scc.SCC_MAX_EXPORT_REDIRECTS + 1)
        self.assertEqual(list(self.folder.iterdir()), [])

    def test_login_redirect_message_distinguishes_it_from_rejected_basic_auth(self):
        redirect = self.response(302)
        redirect.headers["Location"] = "/admin/login/?next=/private"
        with patch.object(scc.requests, "get", return_value=redirect):
            with self.assertRaisesRegex(scc.SCCDownloadError, "login page") as caught:
                scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password")
        self.assertNotIsInstance(caught.exception, scc.SCCAuthenticationError)

    def test_html_login_page_is_not_saved_as_csv(self):
        with patch.object(scc.requests, "get", return_value=self.response(content_type="text/html")):
            with self.assertRaises(scc.SCCDownloadError):
                scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password", "result.csv")
        self.assertEqual(list(self.folder.iterdir()), [])

    def test_interrupted_download_preserves_previous_file_and_hides_error_details(self):
        final = self.folder / "result.csv"
        final.write_bytes(b"old valid data")
        def interrupted():
            yield b"partial"
            raise scc.requests.exceptions.ConnectionError("sensitive exception detail")
        response = self.response()
        response.iter_content.return_value = interrupted()
        with patch.object(scc.requests, "get", return_value=response):
            with self.assertRaises(scc.SCCDownloadError) as caught:
                scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password", final.name)
        self.assertNotIn("sensitive", str(caught.exception))
        self.assertEqual(final.read_bytes(), b"old valid data")
        self.assertEqual(list(self.folder.iterdir()), [final])

    def test_timeout_message_hides_request_exception_details(self):
        with patch.object(scc.requests, "get", side_effect=scc.requests.exceptions.Timeout("sensitive detail")):
            with self.assertRaisesRegex(scc.SCCDownloadError, "timed out") as caught:
                scc.download_with_auth(self.url, str(self.folder), "test-user", "test-password")
        self.assertNotIn("sensitive", str(caught.exception))

    def test_cli_auth_failure_returns_error_without_traceback(self):
        with patch.object(cli, "export_scc_config", side_effect=scc.SCCAuthenticationError("Authentication rejected")):
            with contextlib.redirect_stderr(io.StringIO()) as error:
                result = cli.main(["-c", "713", "-o", str(self.folder)])
        self.assertEqual(result, 1)
        self.assertIn("Authentication rejected", error.getvalue())
        self.assertNotIn("Traceback", error.getvalue())

    def test_cli_cancellation_returns_130(self):
        with patch.object(cli, "export_scc_config", side_effect=KeyboardInterrupt):
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(cli.main(["-c", "713", "-o", str(self.folder)]), 130)


if __name__ == "__main__":
    unittest.main()

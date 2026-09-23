Move these two existing project files into this package without changing their internals:

  utils/parse_intercomparison_args.py -> intercomparison/arguments.py
  utils/intercomparison_stages.py     -> intercomparison/stages.py

The updated __intercomparison_interactive__.py already imports them from these new paths.
No compatibility copies should remain in utils.

import wx

class AtlasUIInputField (wx.Window):
    def __init__ ( self, *args, **kwargs ):
        
        self.label = kwargs.pop ("label", "")
        self.value = kwargs.pop ("value", None)
        self.optional = kwargs.pop ("optional", False)
        self.selected = kwargs.pop ("selected", False)
        self.setting = kwargs.pop ("setting", None)
        self.enabled = kwargs.pop ("enabled", True)
        orient = kwargs.pop("orient", "vertical")
        
        if self.setting is not None:
            value = self.setting.value
            
            if value is None:
                # Setting might be unselected ;)
                value = self.setting.default
                
            self.value = value
            
            self.optional = not self.setting.Mandatory()
            self.selected = self.setting.Selected()
        
        wx.Window.__init__ ( self, *args, **kwargs )
        self.input_window = wx.Window(parent = self)
        
        if orient.lower() == "vertical":
            self.sizer = wx.BoxSizer (wx.VERTICAL)
        else:
            self.sizer = wx.BoxSizer (wx.HORIZONTAL)
        
        if self.optional:
            self.label_ctrl = wx.CheckBox ( parent = self, label = self.label )
            self.label_ctrl.Bind (wx.EVT_CHECKBOX, self.OnCheckBoxClicked)
            self.label_ctrl.SetValue ( self.selected )
        else:
            self.label_ctrl = wx.StaticText ( parent = self, label = self.label )
            
        self.Enable ( self.enabled )
            
        self.sizer.Add ( self.label_ctrl )
        self.sizer.Add ( self.input_window, 1, wx.EXPAND )

        self.SetSizer (self.sizer)
        
        self.update_subscriptions = []
            
    def enable_input_window ( self, enabled = True ):
        # Hack to force re-enabling or re-disabling the window:
        if not self.input_window.Enable (enabled):
            self.input_window.Enable (not enabled)
            self.input_window.Enable (enabled)
        
        if self.setting:
            if enabled:
                self.setting.Select()
            else:
                self.setting.Deselect()
                
    def OnUpdate ( self, func ):
        if not callable (func):
            raise ValueError ("`func` argument is not a callable function!")
            
        if func in self.update_subscriptions:
            raise ValueError ("Method already subscribed to data updates!")
            
        self.update_subscriptions.append ( func )
        
    def RemoveUpdateSubscription ( self, func = None ):
        # if func is None:
            # self.update_subscriptions = []
            
        if func in self.update_subscriptions:
            self.update_subscriptions.remove ( func )
            
    
        
    def OnCheckBoxClicked ( self, evt ):
        self.Select ( evt.IsChecked() )
        
        evt.Skip()
        
    def Enable ( self, enabled = True ):
        self.enabled = enabled
        self.label_ctrl.Enable( enabled )
        self.input_window.Enable ( enabled )
        
        if self.setting:
            self.setting.Select ( enabled and self.selected )
            self.enable_input_window ( enabled )
            
    def Disable ( self ):
        self.Enable ( False )
        
    def Select ( self, selected = True ):
        if not self.optional:
            self.selected = True
            
            if self.setting and not self.enabled:
                self.setting.Select()
            return
            
        self.selected = selected
        
        self.enable_input_window ( self.selected and self.enabled )
        
    def IsSelected ( self ):
        return not self.optional or self.selected
        
    def GetValue ( self ):
        return self.value
        
    def SetValue ( self, value = None ):
        self.value = value
        
    def update_setting ( self, value = None ):
        if value is None:
            value = self.GetValue()
            
        if self.setting:
            self.setting ( value )
            
        for subscription in self.update_subscriptions:
            try:
                subscription( value )
            except Exception:
                continue
        
    @staticmethod
    def FromSetting ( setting, *args, **kwargs ):
        setting_type = type(setting).__name__
        
        if setting_type == 'ATLASSelectOption':
            return CheckBoxField ( setting = setting, *args, **kwargs )
        elif setting_type == 'ATLASTextOption':
            return TextField ( setting = setting, *args, **kwargs )
        elif setting_type == 'ATLASFolderOption':
            return FolderField ( setting = setting, *args, **kwargs )
        elif setting_type == 'ATLASIntegerOption':
            return NumberField ( setting = setting, *args, **kwargs )
        elif setting_type == 'ATLASIntegerRangeOption':
            return RangeField ( setting = setting, *args, **kwargs )
        elif setting_type == 'ATLASDoubleOption':
            return NumberField ( setting = setting, *args, **kwargs )
        elif setting_type == 'ATLASDoubleRangeOption':
            return RangeField ( setting = setting, *args, **kwargs )
        elif setting_type == 'ATLASChoiceOption':
            return SingleChoiceField ( setting = setting, *args, **kwargs )
        elif setting_type == 'ATLASMultipleChoiceOption':
            return MultipleChoiceField ( setting = setting, *args, **kwargs )
            
        return TextField ( setting = setting, *args, **kwargs )
        
class CheckBoxField (AtlasUIInputField):
    def __init__ (self, *args, **kwargs):
        AtlasUIInputField.__init__ ( self, *args, **kwargs )
        
        self.label_ctrl.Hide()
        self.label_ctrl.SetLabel ('')
        
        self.input_sizer = wx.BoxSizer (orient = wx.HORIZONTAL)
        
        self.check = wx.CheckBox ( self.input_window, label = self.label )
        self.Bind ( wx.EVT_CHECKBOX, self.OnCheckChanged, self.check )
        self.input_sizer.Add (self.check, 1, wx.EXPAND)
        
        self.input_window.SetSizer ( self.input_sizer )
        self.SetValue ( self.value )
        self.Select ( True )
        
    def OnCheckChanged ( self, e ):
        self.update_setting()
        e.Skip()
        
    def GetValue ( self ):
        return self.check.IsChecked()
        
    def SetValue (self, value = None):
        if value is None or not isinstance(value, bool):
            self.value = self.check.SetValue (False)
        else:
            self.value = self.check.SetValue(value)
            
        if self.setting:
            self.setting ( self.value )
        
class FolderField (AtlasUIInputField):
    def __init__ ( self, *args, **kwargs ):
        AtlasUIInputField.__init__ ( self, *args, **kwargs )
        
        self.input_sizer = wx.BoxSizer (orient = wx.HORIZONTAL)
        
        self.text = wx.TextCtrl ( self.input_window, style = wx.TE_READONLY )
        
        self.btn = wx.Button ( self.input_window, label = "Browse..." )
        self.btn.Bind ( wx.EVT_BUTTON, self.OnBtnClicked )
        
        self.input_sizer.Add (self.text, 1, wx.EXPAND)
        self.input_sizer.AddSpacer (5)
        self.input_sizer.Add (self.btn, 0, wx.EXPAND)
        
        self.input_window.SetSizer (self.input_sizer)
        
        self.SetValue ( self.value )
        self.Select ( self.IsSelected() )
        
    def OnBtnClicked (self, event):
        dlg = wx.DirDialog ( self, style = wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST )
        if dlg.ShowModal() == wx.ID_OK:
            self.SetValue ( dlg.GetPath() )
        else:
            self.SetValue ( None )
            
        dlg.Destroy()
        
    def GetValue ( self ):
        return self.text.GetValue()
        
    def SetValue (self, value = None):
        final_value = ""
        
        if value is None or not isinstance(value, str):
            self.text.SetValue ("")
            final_value = ""
        else:
            self.text.SetValue(value)
            final_value = value
            
        self.value = final_value
            
        if self.setting:
            self.setting ( final_value )
            
class TextField (AtlasUIInputField):
    def __init__ (self, *args, **kwargs):
        AtlasUIInputField.__init__ ( self, *args, **kwargs )
        
        self.input_sizer = wx.BoxSizer (orient = wx.HORIZONTAL)
        
        self.text = wx.TextCtrl ( self.input_window )
        self.Bind ( wx.EVT_TEXT, self.OnTextChanged, self.text )
        self.input_sizer.Add (self.text, 1, wx.EXPAND)
        
        self.input_window.SetSizer ( self.input_sizer )
        self.SetValue ( self.value )
        self.Select ( self.IsSelected() )
        
    def OnTextChanged ( self, e ):
        self.update_setting()
        e.Skip()
        
    def GetValue ( self ):
        return self.text.GetValue()
        
    def SetValue (self, value = None):
        final_value = ""
        
        if value is None or not isinstance(value, str):
            self.text.SetValue ("")
            final_value = ""
        else:
            self.text.SetValue(value)
            final_value = value
            
        self.value = final_value
            
        if self.setting:
            self.setting ( self.value )
            
class SingleChoiceField (AtlasUIInputField):
    def __init__ (self, *args, **kwargs):
        choices = kwargs.pop ("choices", {})
        AtlasUIInputField.__init__ ( self, *args, **kwargs )
        
        s = kwargs.get("setting", None)
        
        if s and type(s).__name__ == 'ATLASChoiceOption':
            if s.Choices():
                choices = s.Choices()
        
        self.input_sizer = wx.BoxSizer (orient = wx.HORIZONTAL)
        
        self.choice = wx.Choice ( self.input_window, choices = [] )
        self.Bind ( wx.EVT_CHOICE, self.OnChoiceChanged, self.choice )
        self.SetChoices(choices)
        
        self.input_sizer.Add (self.choice, 1, wx.EXPAND)
        
        self.input_window.SetSizer ( self.input_sizer )
        self.SetValue ( self.value )
        self.Select ( self.IsSelected() )
        
    def OnChoiceChanged ( self, event ):
        self.update_setting()
        event.Skip()
        
    def SetChoices ( self, choices = {} ):
        try:
            old_string = self.choice.GetString( self.choice.GetCurrentSelection() )
        except Exception:
            old_string = None
        
        try:
            if choices is None:
                raise ValueError("Choices cannot be None!")
                
            if isinstance(choices, dict):
                self.choices = choices
            elif isinstance(choices, list):
                self.choices = { choice : choice for choice in choices }
                self.choice.Set([choice for choice in self.choices.values()])
            else:
                raise ValueError("Choices not list or dict!")
        except ValueError:
            self.choices = {}

        self.choice.Set([choice for choice in self.choices.values()])
        
        if old_string is not None:
            selection = self.choice.FindString (old_string)
        else:
            selection = wx.NOT_FOUND
        
        if selection == wx.NOT_FOUND:
            self.choice.SetSelection(0)
        else:
            self.choice.SetSelection(selection)
            
    def GetValue ( self ):
        choice = self.choice.GetString ( self.choice.GetSelection() )
        
        for key in self.choices.keys():
            if self.choices[ key ] == choice:
                return key
                
        return None
        
    def SetValue (self, value = None):
        if value is None or not isinstance(value, str):
            self.choice.SetSelection(0)
        else:
            for choice in self.choices.values():
                if self.choices[ value ] == choice:
                
                    if self.setting:
                        self.setting ( value )
                
                    index = self.choice.FindString(choice)
                    if index != wx.NOT_FOUND:
                        self.choice.SetSelection(index)
                    else:
                        self.choice.SetSelection(0)
                        
                    return
                    
        self.choice.SetSelection(0)
        
class MultipleChoiceDialog(wx.Dialog):
    def __init__ (self, parent, title="Select options", choices = {}, selected = []):
        super(MultipleChoiceDialog, self).__init__(None, title = title)
        
        self.choices = choices
        self.initial_selection = selected
        self.selected = selected
        
        self.SetSize((-1, 400))
        
        panel = wx.Panel(self)
        panelVBox=wx.BoxSizer(wx.VERTICAL)
        
        TopHBox=wx.BoxSizer(wx.HORIZONTAL)
        
        self.selectAllBtn = wx.Button(panel, label="Select all")
        self.selectAllBtn.Bind(wx.EVT_BUTTON, self.OnSelectAllBtnClicked)
        self.selectNoneBtn = wx.Button(panel, label="Select none")
        self.selectNoneBtn.Bind(wx.EVT_BUTTON, self.OnSelectNoneBtnClicked)
        
        TopHBox.Add(self.selectAllBtn, proportion=0, flag=wx.RIGHT, border=5)
        TopHBox.Add(self.selectNoneBtn, proportion=0, flag=wx.LEFT, border=5)
        
        self.list = wx.ListBox(panel, style=wx.LB_EXTENDED)
        
        panelVBox.Add(TopHBox, proportion=0, flag=wx.ALL, border=10)
        panelVBox.Add(self.list, proportion=1, flag=wx.EXPAND)
        
        panel.SetSizer(panelVBox)
        
        BottomHBox=wx.BoxSizer(wx.HORIZONTAL)
        
        self.applyBtn = wx.Button(self, label="OK")
        self.applyBtn.Bind(wx.EVT_BUTTON, self.OnApplyBtnClicked)
        self.cancelBtn = wx.Button(self, label="Cancel")
        self.cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancelBtnClicked)
        
        BottomHBox.Add(self.applyBtn, proportion=0, flag=wx.RIGHT, border=5)
        BottomHBox.Add(self.cancelBtn, proportion=0, flag=wx.LEFT, border=5)
        
        VBox=wx.BoxSizer(wx.VERTICAL)
        VBox.Add(panel, proportion=1, flag=wx.EXPAND|wx.BOTTOM, border=10)
        VBox.Add(BottomHBox, flag=wx.CENTER|wx.ALL, border=10)
        
        self.SetSizer(VBox)
        
        self.PopulateListBox ()
        
    def OnSelectAllBtnClicked(self, e):
        for i in range(self.list.GetCount()):
            self.list.SetSelection(i)
            
    def OnSelectNoneBtnClicked(self, e):
        for i in range(self.list.GetCount()):
            self.list.Deselect(i)
            
    def OnApplyBtnClicked(self, e):
        names=self.list.GetItems()
        selected_names = []
        
        self.selected = []
        
        for i in range(self.list.GetCount()):
            checked=self.list.IsSelected(i)
            name = names[i]
            
            if checked:
                selected_names.append ( name )
                
        for name in selected_names:
            for choice in self.choices.keys():
                if self.choices[ choice ] == name:
                    self.selected.append (choice)
        
        self.EndModal(wx.ID_OK)
            
    def OnCancelBtnClicked(self, e):
        self.selected = self.initial_selection
        self.EndModal(wx.ID_CANCEL)
        
    def PopulateListBox (self):
        # Insert items:
        items = [ text for text in self.choices.values() ]
        
        if len(items):
            self.list.InsertItems ( items, 0 )
            self.list.EnsureVisible(0)
        
        # Also select appropiate items:
        for selected in self.selected:
            display_name = self.choices.get(selected, None)
            
            if display_name is not None:
                self.list.SetStringSelection(display_name, True)
                
    def GetSelectedChoices(self):
        return self.selected
        
        
class MultipleChoiceField (AtlasUIInputField):
    def __init__ ( self, *args, **kwargs ):
        choices = kwargs.pop ("choices", {})
        
        s = kwargs.get("setting", None)
        if s and type(s).__name__ == 'ATLASMultipleChoiceOption':
            if s.Choices():
                choices = s.Choices()
            
        self.SetChoices(choices)
        
        AtlasUIInputField.__init__ ( self, *args, **kwargs )
        
        self.input_sizer = wx.BoxSizer (orient = wx.HORIZONTAL)
        
        self.text = wx.TextCtrl ( self.input_window, style = wx.TE_READONLY )
        
        self.btn = wx.Button ( self.input_window, label = "Select..." )
        self.btn.Bind ( wx.EVT_BUTTON, self.OnBtnClicked )
        
        self.input_sizer.Add (self.text, 1, wx.EXPAND)
        self.input_sizer.AddSpacer (5)
        self.input_sizer.Add (self.btn, 0, wx.EXPAND)
        
        self.input_window.SetSizer (self.input_sizer)
        
        self.SetValue ( self.value )
        self.Select ( self.IsSelected() )
        
    def OnBtnClicked (self, event):
        dlg = MultipleChoiceDialog (self, title=self.label, choices = self.choices, selected = self.selected_choices)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.selected_choices = dlg.GetSelectedChoices()
            
        dlg.Destroy()
        # Refresh text field!
        self.SetValue ([selected for selected in self.selected_choices])
        
        if self.setting:
            self.setting ( self.selected_choices )
            
    def SetChoices ( self, choices = {} ):
        if choices is None or not isinstance(choices, dict):
            self.choices = {}
        else:
            self.choices = choices
            
        self.selected_choices = []
        
    def GetValue ( self ):
        return self.selected_choices
        
    def SetValue (self, value = None):
        self.selected_choices = []
        
        if value is not None and isinstance(value, list):
            for selected_value in value:
                if selected_value in self.choices.keys():
                    self.selected_choices.append( selected_value )
                    
        if self.setting:
            self.setting ( self.selected_choices )
                    
        selected_choices_labels = [ self.choices[choice] for choice in self.selected_choices ]
        
        if len(selected_choices_labels):
            self.text.SetValue(", ".join(selected_choices_labels))
        else:
            self.text.SetValue("None selected")
            
class NumberField (AtlasUIInputField):
    def __init__ (self, *args, **kwargs):
        self.min = kwargs.pop("min", 0)
        self.max = kwargs.pop("max", 100)
        self.increment = kwargs.pop("increment", 1)
        
        s = kwargs.get ("setting", None)
        
        if s and (type(s).__name__ == 'ATLASDoubleOption' or type(s).__name__ == 'ATLASIntegerOption'):
            self.min = s.Min()
            self.max = s.Max()
            self.increment = s.Inc()
        
        AtlasUIInputField.__init__ ( self, *args, **kwargs )
        
        self.input_sizer = wx.BoxSizer (orient = wx.HORIZONTAL)
        
        self.spinner = wx.SpinCtrlDouble ( self.input_window )
        
        if self.min:
            self.spinner.SetMin ( self.min )
            
        if self.max:
            self.spinner.SetMax ( self.max )
            
        if self.increment:
            self.spinner.SetIncrement ( self.increment )
        
        self.Bind ( wx.EVT_SPINCTRLDOUBLE, self.OnValueChanged, self.spinner )
        self.input_sizer.Add (self.spinner, 1, wx.EXPAND)
        
        self.input_window.SetSizer ( self.input_sizer )
        self.SetValue ( self.value )
        self.Select ( self.IsSelected() )
        
    def OnValueChanged ( self, event ):
        value = self.spinner.GetValue()
        
        if self.increment == 1:
            value = int (value)
            
        self.update_setting ( value )
            
        event.Skip()
        
    def GetValue ( self ):
        return self.spinner.GetValue()
        
    def SetValue (self, value = None):
        if value:
            self.spinner.SetValue (value)
            
    def SetIncrement ( self, increment ):
        self.increment = increment
        self.spinner.SetIncrement ( self.increment )
            
class RangeField (AtlasUIInputField):
    def __init__ (self, *args, **kwargs):
        self.min = kwargs.pop("min", None)
        self.max = kwargs.pop("max", None)
        self.increment = kwargs.pop("increment", 1)
        
        s = kwargs.get ("setting", None)
        
        if s and (type(s).__name__ == 'ATLASDoubleRangeOption' or type(s).__name__ == 'ATLASIntegerRangeOption'):
            self.min = s.Min()
            self.max = s.Max()
            self.increment = s.Inc()
        
        AtlasUIInputField.__init__ ( self, *args, **kwargs )
        
        self.input_sizer = wx.BoxSizer (orient = wx.HORIZONTAL)
        
        self.start_spinner = wx.SpinCtrlDouble ( self.input_window )
        self.end_spinner = wx.SpinCtrlDouble ( self.input_window )
        
        if self.min:
            self.start_spinner.SetMin ( self.min )
            
        if self.max:
            self.end_spinner.SetMax ( self.max )
            
        if self.increment:
            self.start_spinner.SetIncrement ( self.increment )
            self.end_spinner.SetIncrement ( self.increment )
        
        dash_label = wx.StaticText ( self.input_window, label = "-" )
        
        self.Bind ( wx.EVT_SPINCTRLDOUBLE, self.OnValueChanged, self.start_spinner )
        self.Bind ( wx.EVT_SPINCTRLDOUBLE, self.OnValueChanged, self.end_spinner )
        
        self.input_sizer.Add (self.start_spinner, 1, wx.EXPAND)
        self.input_sizer.AddSpacer (5)
        self.input_sizer.Add (dash_label)
        self.input_sizer.AddSpacer (5)
        self.input_sizer.Add (self.end_spinner, 1, wx.EXPAND)
        
        self.input_window.SetSizer ( self.input_sizer )
        self.SetValue ( self.value )
        self.Select ( self.IsSelected() )
        
    def OnValueChanged ( self, event ):
        if self.min is not None:
            self.start_spinner.SetRange ( self.min, self.end_spinner.GetValue() )
            
        if self.max is not None:
            self.end_spinner.SetRange ( self.start_spinner.GetValue(), self.max )
            
        self.update_setting()
        event.Skip()
        
    def GetValue ( self ):
        if self.increment == 1:
            value = (int(self.start_spinner.GetValue()), int(self.end_spinner.GetValue()))
        else:
            value = (self.start_spinner.GetValue(), self.end_spinner.GetValue())
            
        return value
        
    def SetValue (self, value = None):
        if value is None or not isinstance(value, tuple):
            try:
                self.start_spinner.SetValue (self.min)
                self.end_spinner.SetValue (self.max)
            except Exception:
                pass
        elif self.min and self.max:
            start_value = min ( max ( self.min, value[0] ), self.max )
            end_value = min ( max ( self.min, value[1] ), self.max )
            
            self.start_spinner.SetValue( start_value )
            self.end_spinner.SetValue( end_value )
        else:
            self.start_spinner.SetValue ( value[0] )
            self.end_spinner.SetValue ( value[1] )
            
    def SetIncrement ( self, increment ):
        self.increment = increment
        self.start_spinner.SetIncrement ( self.increment )
        self.end_spinner.SetIncrement ( self.increment )
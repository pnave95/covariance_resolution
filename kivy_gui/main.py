from kivy.uix.tabbedpanel import TabbedPanel
from kivy.app import App
from kivy.uix.floatlayout import FloatLayout
from kivy.uix.tabbedpanel import TabbedPanel
from kivy.uix.tabbedpanel import TabbedPanelHeader
from kivy.uix.label import Label
from kivy.uix.button import Button
from kivy.uix.gridlayout import GridLayout

from kivy.uix.textinput import TextInput

from kivy.properties import ObjectProperty

from kivy.lang import Builder

Builder.load_string("""

<LatticeInput>:
    size_hint: 1, 1
    pos_hint: {'center_x': 0.5, 'center_y': 0.5}

    a11: a11
    a12: a12
    a13: a13
    a21: a21
    a22: a22
    a23: a23
    a31: a31
    a32: a32
    a33: a33
    GridLayout:
        rows: 3
        cols: 3

        padding: 5
        spacing: 5

        TextInput:
            id: a11
        TextInput: 
            id: a12
        TextInput:
            id: a13
        TextInput:
            id: a21
        TextInput:
            id: a22
        TextInput:
            id: a23
        TextInput:
            id: a31
        TextInput:
            id: a32
        TextInput:
            id: a33

        


<RootTabPanel>:
    size_hint: 1, 1
    pos_hint: {'center_x': .5, 'center_y': .5}
    do_default_tab: False
    tab_width: 150

    TabbedPanelItem:
        text: 'Lattice Parameters'
        LatticeInput:
        

    TabbedPanelItem:
        text: 'Instrument Settings'
        GridLayout:
            rows: 2
            cols: 2

            padding: 10
            spacing: 10
            Label:
                text: 'Instrument Settings'
            Button:
                text: 'Button that does nothing'
            Label:
                text: '(More instrument Settings)'
            Button:
                text: 'Another useless button'
    TabbedPanelItem:
        text: 'Resolution Analysis'
        RstDocument:
            text:
                '\\n'.join(("Hello world", "-----------",
                "You are in the third tab."))


""")

class LatticeInput(GridLayout):
    a11 = ObjectProperty(None)
    a12 = ObjectProperty(None)
    a13 = ObjectProperty(None)
    a21 = ObjectProperty(None)
    a22 = ObjectProperty(None)
    a23 = ObjectProperty(None)
    a31 = ObjectProperty(None)
    a32 = ObjectProperty(None)
    a33 = ObjectProperty(None)
    #pass

class RootTabPanel(TabbedPanel):
    #txt_inpt = ObjectProperty(None)
    pass


class TabbedPanelApp(App):
    def build(self):
        return RootTabPanel()

#class RootWidget(FloatLayout):
    #pass



# class TabbedPanelApp(App):
#       def build(self):
#           tb_panel= TabbedPanel()
 
#           #Create text tab          
#           th_text_head = TabbedPanelHeader(text='Lattice Parameters')
#           th_text_head.content= Label(text='This is my text content')
 
#           #Create image tab
#           #th_img_head= TabbedPanelHeader(text='Image tab')
#           #th_img_head.content= Image(source='sample.jpg',pos=(400, 100), size=(400, 400))
 
#           #Create button tab
#           th_btn_head = TabbedPanelHeader(text='Button tab')
#           th_btn_head.content= Button(text='This is my button',font_size=20)
 
#           tb_panel.add_widget(th_text_head)
#           #tb_panel.add_widget(th_img_head)
#           tb_panel.add_widget(th_btn_head)          
 
#           return tb_panel

#class MainApp(App):
#    def build(self):
#        return RootWidget()

if __name__ == '__main__':
    #MainApp().run()
    TabbedPanelApp().run()

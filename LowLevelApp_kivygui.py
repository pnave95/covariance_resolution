from kivy.app import App
from kivy.uix.widget import Widget
from kivy.uix.gridlayout import GridLayout
from kivy.uix.label import Label


class InputScreen(GridLayout):

	def __init__(self, **kwargs):
		super(InputScreen, self).__init__(**kwargs)
		self.cols = 2
		#self.add_widget(Label(text=''))


class LowLevelApp(App):

	def build(self):
		return InputScreen()


if __name__ == '__main__':
	LowLevelApp.run()
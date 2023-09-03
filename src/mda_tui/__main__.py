from textual import app, widgets


class MDA(app.App[str]):
    def __init__(self):
        super().__init__()
        self.title = "MDA TUI"
        self.sub_title = "On-the-fly transformations and trajectory analysis in your terminal!"

    def compose(self) -> app.ComposeResult:
        yield widgets.Header()
        yield widgets.Footer()


def main():
    app = MDA()
    app.run()

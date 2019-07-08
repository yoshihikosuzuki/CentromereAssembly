from flask import Flask
from .read_viewer import ReadViewer


@app.route("/")
def index():
    # Select db_file, las_file, out_dir, gepard via Window
    # TODO: save a viewer setting to a file?

    # Create a viewer
    v = ReadViewer()

    # Select a read

    # Show plot
    v.show()

    return

import dash
from dash import dcc
from dash import html
from dash import callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State

# To open in browser
import webbrowser
from threading import Timer

# Upload management
import os
import base64
import datetime
import io
from pathlib import Path

# Local Imports
from slab_gen import Slab, InterfaceBuilder
from .cell_match_strain import gen_plot

class App:
    def __init__(self):
        # Goal to insert desired fig generated from chosen files
        self.app = dash.Dash(
            "Interactive Interface Graphing",
            external_stylesheets=[dbc.themes.BOOTSTRAP],
        )
        self.app.layout = self.build_startup_html()

        # Upload Management
        self.UPLOAD_DIRECTORY = self.setup_upload_dir()

        # Startup
        Timer(1, self.open_in_browser).start()
        self.app.run_server(debug=True, port=7000, use_reloader=False)

    # Building starttup page where you can select graph to genera
    def build_startup_html(self):
        return html.Div(
            [
                html.H1("Gen_Plot Interactive Plot - Adelstein Research Group"),
                html.Small("By Arye Oskotsky"),
                html.H3("Select structures to plot"),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    [
                                        dcc.Upload(
                                            id="slab1_upload",
                                            children=html.Button(
                                                "Upload Structure - Rotated"
                                            ),
                                            multiple=False,
                                        )
                                    ],
                                    className="col",
                                ),
                                html.Div(
                                    [
                                        dcc.Upload(
                                            id="slab2_upload",
                                            children=html.Button(
                                                "Upload Structure - Static"
                                            ),
                                            multiple=False,
                                        )
                                    ],
                                    className="col",
                                ),
                            ],
                            className="row",
                        )
                    ],
                    className="container",
                ),
                html.Div([dcc.Graph('graph_content')], id='graph_container')
            ]
        )
    
    # Generate the graphs. Current iteration will generate all graphs. 
    # Currently:
    #          1. # Atoms & Energy & Strain
    @callback(
        Output('graph_container', 'children'),
        [
            Input('slab1_upload', 'contents'),
            Input('slab1_upload', 'filename'),
            Input('slab2_upload', 'contents'),
            Input('slab2_upload', 'filename')
        ],
        State('gen_plot', 'value'),
        prevent_intial_call=True
    )
    def generate_plots(self, slab1_upload, slab1_fname, slab2_upload, slab2_fname, plot_selector_dropdown):
        if (slab1_upload is None or slab2_upload is None): return dash.no_update
        slab1 = self.save_file(slab1_upload, slab1_fname)
        slab2 = self.save_file(slab2_upload, slab1_fname)
        interface = InterfaceBuilder(slab1, slab2)

        return gen_plot(interface)


        

    # Force open in browser
    def open_in_browser(self):
        if not os.environ.get("WERKZUEG_RUN_MAIN"):
            webbrowser.open("http://localhost:{}".format(7000))

    # Saves file and returns path to file
    def save_file(self, file : str, filename : str) -> str:
        _, content_string = file.split(',')
        decoded = base64.b64decode(content_string)
        fpath = os.path.join(self.UPLOAD_DIRECTORY, filename)
        
        with open(fpath, 'wb') as f:
            f.write(decoded)
        
        return fpath

    # Sets up upload file if it is not already set up and saves abs path to uploads
    def setup_upload_dir(self):
        if not os.path.exists('uploads'):
            os.mkdir('uploads')
        
        return Path('uploads').resolve()

if __name__ == "__main__":
    App()

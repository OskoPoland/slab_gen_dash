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
from utils import energy_nat_strain

class App:
    def __init__(self):
        # Interface Data Variables
        self.slab1 = None
        self.slab2 = None
        self.interface = None
        
        # Goal to insert desired fig generated from chosen files
        self.app = dash.Dash(
            "Interactive Interface Graphing",
            external_stylesheets=[dbc.themes.BOOTSTRAP],
        )
        self.app.layout = self.build_startup_html()

        # Upload Management
        self.UPLOAD_DIRECTORY = self.setup_upload_dir()
    
        if self.app is not None and hasattr(self, "callbacks"):
            self.callbacks(self.app)

        # Startup
        self.start()
    
    # Scoping callback methods as static methods within class method
    def callbacks(self, _app):
        @_app.callback(
        [
            Output('slab1_upload', 'title'),
            Output('slab2_upload', 'title')
        ],
        [
            Input('slab1_upload', 'contents'),
            Input('slab1_upload', 'filename'),
            Input('slab2_upload', 'contents'),
            Input('slab2_upload', 'filename'),
        ]
        )
        def upload_structures(slab1_upload, slab1_fname, slab2_upload, slab2_fname):
            if (slab1_upload is None or slab2_upload is None): 
                return dash.no_update
                
            self.slab1 = self.save_file(slab1_upload, slab1_fname)
            self.slab2 = self.save_file(slab2_upload, slab1_fname)
            self.interface = InterfaceBuilder(self.slab1, self.slab2)
            
            static_structure = slab1_fname.split("\\")[-1]
            nonstatic_structures = slab2_fname.split("\\")[-1]
            
            return static_structure, nonstatic_structures
        
        
        @_app.callback(
            # Output('graph_content', 'figure'),
            Output('graph_content', 'figure'),
            [Input('gen_graphs', 'n_clicks')]
        )
        def generate_graphs(n_clicks):
            if self.interface is None:
                return dash.no_update
            
            # Isolate id's of changed properties since last callback trigger
            changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
            
            if 'gen_graphs' in changed_id:
                # Button was pressed
                return energy_nat_strain(self.interface)
            else:
                return dash.no_update
    
    def start(self):
        print("Starting Dash Application")
        Timer(1, self.open_in_browser).start()
        self.app.run_server(debug=True, port=7000, use_reloader=False)
        print("App Has Started")

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
                                        ),
                                        html.Div(id='click_val')
                                    ],
                                    className="col",
                                ),
                                html.Div(
                                    [
                                        html.Button(
                                            "Generate Graphs",
                                            id="gen_graphs"
                                        )
                                    ]
                                )
                            ],
                            className="row",
                        )
                    ],
                    className="container",
                ),
                html.Div([dcc.Graph('graph_content')], id='graph_container')
            ]
        )
    
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

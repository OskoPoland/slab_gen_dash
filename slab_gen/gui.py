from tkinter import *
import tkinter as tk
from tkinter.filedialog import askopenfilename, askdirectory
from .gen_interface import parse_args
import os, yaml

#TODO: Creating interactive plotly plots
class App(tk.Frame):
    def __init__(self, master, root):
        super().__init__(master)
        # self.pack()
        self.master = master
        self.leftparameterframe = tk.Frame(root)
        self.rightparameterframe = tk.Frame(root)
        self.submitframe = tk.Frame(root)
        self.entryfields = {}
        self.fieldvars = {}
                        
        #-- Generate fields
        self.fieldvars.update({'structure1' : tk.StringVar()})
        self.entryfields.update({
            'structure1' : [
                Label(master, text='Path to Structure1').grid(row=1,column=1), 
                Entry(master, width=80, textvariable=self.fieldvars.get('structure1')),
                Button(master, text='Choose', command=lambda: self.select_location('structure1')).grid(row=1,column=3),
            ]
        })
        
        self.fieldvars.update({'structure2' : tk.StringVar()})
        self.entryfields.update({
            'structure2' : [
                Label(master, text='Path to Structure2').grid(row=2,column=1),
                Entry(master, width=80, textvariable=self.fieldvars.get('structure2')),
                Button(master, text='Choose', command=lambda: self.select_location('structure2')).grid(row=2,column=3),
            ]
        })

        self.fieldvars.update({'output' : tk.StringVar()})
        self.entryfields.update({
            'output' : [
                Label(master, text='Output Directory').grid(row=3,column=1),
                Entry(master, width=80, textvariable=self.fieldvars.get('output')),
                Button(master, text='Choose', command=lambda: self.select_location('output')).grid(row=3,column=3),
            ]
        })

        self.fieldvars.update({'config' : tk.StringVar()})
        self.entryfields.update({
            'config' : [
                Label(master, text='Config Location').grid(row=4,column=1),
                Entry(master, width=80, textvariable=self.fieldvars.get('config')),
                Button(master, text='Choose', command=lambda: self.select_location('config')).grid(row=4,column=3)
            ]
        })
        
        self.fieldvars.update({'ncand' : tk.StringVar()})
        self.entryfields.update({
            'ncand' : [
                Label(self.leftparameterframe, text='Number of Candidates').grid(row=4,column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('ncand')),
            ]
        })
        
        self.fieldvars.update({'hkl1' : tk.StringVar()})
        self.entryfields.update({
            'hkl1' : [
                Label(self.leftparameterframe, text='Structure 1 Miller Indicies (Format: X Y Z)').grid(row=5, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('hkl1')),
            ]
        })
        
        self.fieldvars.update({'hkl2' : tk.StringVar()})
        self.entryfields.update({
            'hkl2' : [
                Label(self.leftparameterframe, text='Structure 2 Miller Indicies (Format: X Y Z)').grid(row=6, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('hkl2')),
            ]
        })
        
        self.fieldvars.update({'nlay' : tk.StringVar()})
        self.entryfields.update({
            'nlay' : [
                Label(self.leftparameterframe, text='Number of Layers').grid(row=7, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('nlay')),
            ]
        })
        
        self.fieldvars.update({'vac_height' : tk.StringVar()})
        self.entryfields.update({
            'vac_height' : [
                Label(self.leftparameterframe, text='Vaccuum Height').grid(row=8, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('vac_height')),
            ]
        })
        
        self.fieldvars.update({'nmax' : tk.StringVar()})
        self.entryfields.update({
            'nmax' : [
                Label(self.leftparameterframe, text='NMAX').grid(row=9, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('nmax')),
            ]
        })

        self.fieldvars.update({'mmax' : tk.StringVar()})
        self.entryfields.update({
            'mmax' : [
                Label(self.leftparameterframe, text='MMAX').grid(row=10, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('mmax')),
            ]
        })
        
        self.fieldvars.update({'theta_min' : tk.StringVar()})
        self.entryfields.update({
            'theta_min' : [
                Label(self.leftparameterframe, text='Theta Min').grid(row=11, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('theta_min')),
            ]
        })
        
        self.fieldvars.update({'theta_max' : tk.StringVar()})
        self.entryfields.update({
            'theta_max' : [
                Label(self.leftparameterframe, text='Theta Max').grid(row=12, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('theta_max')),
            ]
        })
        
        self.fieldvars.update({'dtheta' : tk.StringVar()})
        self.entryfields.update({
            'dtheta' : [
                Label(self.leftparameterframe, text='D-Theta').grid(row=13, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('dtheta')),
            ]
        })
        
        self.fieldvars.update({'strain_cutoff' : tk.StringVar()})
        self.entryfields.update({
            'strain_cutoff' : [
                Label(self.leftparameterframe, text='Strain Cutoff').grid(row=14, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('strain_cutoff')),
            ]
        })
        
        self.fieldvars.update({'theta_cutoff' : tk.StringVar()})
        self.entryfields.update({
            'theta_cutoff' : [
                Label(self.leftparameterframe, text='Theta Cutoff').grid(row=15, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('theta_cutoff')),
            ]
        })
        
        self.fieldvars.update({'strain_weight' : tk.StringVar()})
        self.entryfields.update({
            'strain_weight' : [
                Label(self.leftparameterframe, text='Strain Weight').grid(row=16, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('strain_weight')),
            ]
        })
        
        self.fieldvars.update({'area_weight' : tk.StringVar()})
        self.entryfields.update({
            'area_weight' : [
                Label(self.leftparameterframe, text='Area Weight').grid(row=17, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('area_weight')),
            ]
        })

        self.fieldvars.update({'angle_weight' : tk.StringVar()})
        self.entryfields.update({
            'angle_weight' : [
                Label(self.leftparameterframe, text='Angle Weight').grid(row=18, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('angle_weight')),
            ]
        })
        
        self.fieldvars.update({'grid_spacing' : tk.StringVar()})
        self.entryfields.update({
            'grid_spacing' : [
                Label(self.leftparameterframe, text='Grid Spacing').grid(row=19, column=1,sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('grid_spacing')),
            ]
        })

        self.fieldvars.update({'zshift' : tk.StringVar()})
        self.entryfields.update({
            'zshift' : [
                Label(self.leftparameterframe, text='Z-Shift').grid(row=20, column=1, sticky='w'),
                Entry(self.leftparameterframe, width=20, textvariable=self.fieldvars.get('zshift'))
            ]
        })
        
        self.fieldvars.update({'save_config' : tk.BooleanVar()})
        self.entryfields.update({
            'save_config' : [
                Checkbutton(self.rightparameterframe,
                            text="Save Config",
                            onvalue=True,
                            variable=self.fieldvars.get('save_config'),
                            offvalue=False).grid(row=2,column=1,sticky='w'),
            ]
        })
        
        self.fieldvars.update({'view' : tk.BooleanVar()})
        self.entryfields.update({
            'view' : [
                Checkbutton(self.rightparameterframe,
                            text="View Unit Cells",
                            onvalue=True,
                            variable=self.fieldvars.get('view'),
                            offvalue=False).grid(row=3,column=1,sticky='w'),     
            ]
        })

        self.fieldvars.update({'gen_trans' : tk.BooleanVar()})          
        self.entryfields.update({
            'gen_trans' : [
                Checkbutton(self.rightparameterframe,
                            text="Generate Transitions",
                            onvalue=True,
                            variable=self.fieldvars.get('gen_trans'),
                            offvalue=False).grid(row=4,column=1,sticky='w'),     
            ]
        })
        
        self.fieldvars.update({'verbose' : tk.BooleanVar()})
        self.entryfields.update({
            'verbose' : [
                Checkbutton(self.rightparameterframe,
                            text="Verbose",
                            onvalue=True,
                            offvalue=False,
                            variable=self.fieldvars.get('verbose'))
                            .grid(row=6,column=1,sticky='w')
            ]
        })
        
        self.entryfields.update({
            'submit' : [
                Button(self.submitframe,
                       text='Submit', 
                       command=lambda: self.run_gen_interface()).grid(row=1, column=1)
            ]
        })

        self.textout = Text(self.rightparameterframe, 
                            width=50, 
                            height=15,
                            wrap='word',
                            state='disabled')
        self.textout.tag_configure('err', foreground='red')
        self.textout.tag_configure('def', foreground='black')
        
        self.textclear = Button(self.submitframe, text='clear', command=lambda: self.clear_err())
        
        

        self.pack_all()
        
#TODO: Test config population 
#TODO: Errors on checking boolean typed values
    def select_location(self, bound_button: str) -> None:
        """Populates entry fields upon file/directory selection"""
        if (bound_button == 'output'):
            file_path = askdirectory(mustexist=True)
        else: 
            file_path = askopenfilename()
            
        bound_field = self.entryfields.get(bound_button)[1]
        bound_field.delete(0, tk.END)
        bound_field.insert(0, file_path)
        bound_field.config(fg='black')

        #If config was populated then known fields need to be populated
        if bound_button == 'config':
            config_path = self.entryfields.get(bound_button)[1].get()

            #Another check that config exists and can be safe loaded
            if os.path.exists(config_path):
                with open(config_path, 'r') as cfg:
                    config = yaml.safe_load(cfg)

                    if config is None:
                        self.write_criterr("Config path is valid but could not be loaded. Check that it is a .yaml file.")
                        return

                    #Iterate over keys in config and use the matching value for matching text variables
                    for key in config.keys():
                        key_val = config.get(key)
                        field_to_pop = self.fieldvars.get(key)

                        #Requires string conversion of string bools to bools
                        if key_val == "True":
                            field_to_pop.set(True)
                        elif key_val == "False": 
                            field_to_pop.set(False)
                        else:
                            field_to_pop.set(key_val)
            else:
                self.write_warn('Invalid Config Path: ' + config_path + '. Empty cells will use default values.')
    
    def run_gen_interface(self) -> None:
        """
        Take populated fields and generate an argument string for gen_interface.py. Emtpy parameter fields will 
        be populated with default values in the argument string.
        """
        final_args = []
        path_args = ['output', 'structure1', 'structure2', 'config']
        flag_args = ['view', 'gen_trans', 'save_config', 'verbose']
        
        for args in self.fieldvars.keys():
            if args in path_args:
                fpath = self.fieldvars.get(args).get()
                if os.path.exists(fpath):
                    final_args.extend(['--' + args, fpath])
                elif args == 'output':
                    outpath = self.fieldvars.get('output').get()
                    if not os.path.exists(outpath):
                        self.write_criterr('Defaulting to CWD - ' + os.getcwd())
                        final_args.extend(['-o', os.getcwd()])
                    else:
                        final_args.extend(['-o', args.output])
                elif args == 'config':
                    config_path = self.fieldvars.get('config').get()
                    if not os.path.exists(config_path):
                        self.write_criterr('Config path does not exist! Cannot load values')
                else:
                    self.write_criterr('Destination ' + args + ' does not exist.')
            elif args in flag_args:
                if self.fieldvars.get(args).get():
                    final_args.append('--' + args)
            else:
                value = self.fieldvars.get(args).get()
                if args == 'hkl1' or args == 'hkl2' and value != "":
                    indices = value.split(' ')
                    if len(indices) == 3:
                        final_args.extend(['--' + args, indices[0], indices[1], indices[2]])
                    else:
                        self.write_criterr('Incorrect Miller Inidice Format')
                #Last is trick to determine if string is formated like a float
                elif value != "" and (value.isnumeric() or value.replace('.', '').isnumeric()):
                    final_args.extend(['--' + args, value])
                else:
                    self.write(args + ' using default values')
        
        parse_args(final_args)
    
    def write_warn(self, towrite):
        """Writes warnings to text box"""
        self.textout.config(state='normal')
        self.textout.insert(tk.END, 'WARNING: ' + towrite + '\n', 'def')
        self.textout.config(state='disabled')

    def write_criterr(self, towrite):
        """Writes errors to text box"""
        self.textout.config(state='normal')
        self.textout.insert(tk.END, 'CRITICAL ERROR: ' + towrite + '\n', 'err')
        self.textout.config(state='disabled')

    def write(self, towrite):
        """Writes out information to text box"""
        self.textout.config(state='normal')
        self.textout.insert(tk.END, 'INFO: ' + towrite + '\n')
        self.textout.config(state='disabled')

    def clear_err(self) -> None:
        """Clear error out text box"""
        self.textout.config(state='normal')
        self.textout.delete('1.0', tk.END)
        self.textout.config(state='disabled')

        for args in self.fieldvars:
            self.fieldvars.get(args).set("")
        
    def pack_all(self) -> None:
        """Packs all text entry fields and frames"""
        self.entryfields.get('structure1')[1].grid(row=1, column=2)
        self.entryfields.get('structure2')[1].grid(row=2, column=2)
        self.entryfields.get('output')[1].grid(row=3, column=2)
        self.entryfields.get('config')[1].grid(row=4, column=2)
        
        #Place parameters
        self.leftparameterframe.grid(row=5,column=1,pady=5)
        self.rightparameterframe.grid(row=5,column=2,pady=5,sticky='nw')
        self.submitframe.grid(row=6, column=1)
        
        self.entryfields.get('ncand')[1].grid(row=4, column=2)
        self.entryfields.get('hkl1')[1].grid(row=5, column=2)
        self.entryfields.get('hkl2')[1].grid(row=6, column=2)
        self.entryfields.get('nlay')[1].grid(row=7, column=2)
        self.entryfields.get('vac_height')[1].grid(row=8, column=2)
        self.entryfields.get('nmax')[1].grid(row=9, column=2)
        self.entryfields.get('mmax')[1].grid(row=10, column=2)
        self.entryfields.get('theta_min')[1].grid(row=11, column=2)
        self.entryfields.get('theta_max')[1].grid(row=12, column=2)
        self.entryfields.get('dtheta')[1].grid(row=13, column=2)
        self.entryfields.get('strain_cutoff')[1].grid(row=14, column=2)
        self.entryfields.get('theta_cutoff')[1].grid(row=15, column=2)
        self.entryfields.get('strain_weight')[1].grid(row=16, column=2)
        self.entryfields.get('area_weight')[1].grid(row=17, column=2)
        self.entryfields.get('angle_weight')[1].grid(row=18, column=2)
        self.entryfields.get('grid_spacing')[1].grid(row=19, column=2)
        self.entryfields.get('zshift')[1].grid(row=20, column=2)

        self.textout.grid(row=7, column=1, padx=10)


def run():
    root = tk.Tk()
    
    app = App(root, root)
    app.master.title("Interface Generator")
    app.master.minsize(500,500)
    app.mainloop()


        

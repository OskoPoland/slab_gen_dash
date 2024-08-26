import plotly.express as px
import pandas as pd
import numpy as np
from slab_gen import InterfaceBuilder, Slab

test_chem1 = "C:\\Users\\Vinr5\\OneDrive\\Documents\\Repositories\\Adelstein-research\\codebase\\slab_gen-Arye\\scripts\\Build_interface\\Structures\\Li2CO3.POSCAR.vasp"
test_chem2 = "C:\\Users\\Vinr5\\OneDrive\\Documents\\Repositories\\Adelstein-research\\codebase\\slab_gen-Arye\\scripts\\Build_interface\\Structures\\LiF.POSCAR.vasp"

def gen_plot():
    slab1 = Slab(test_chem1, [1,0,0])
    slab2 = Slab(test_chem2, [1,0,0])
    ib = InterfaceBuilder(slab1, slab2)
    
    # y = number of atoms - nat
    # x = mean abs strain - eps_av
    # Color = log_10(strain scores) - score
    
    # candidates contains nested dictionaries. To convert to df must be flattened
    candidates = ib.get_all_candidates()
    for cand in candidates:
        eps = cand['eps']
        del cand['eps']
        cand['eps_ave'] = eps['eps_ave']
        
    candidates_df = pd.DataFrame.from_dict(ib.get_all_candidates())

    # Build plotly scatter plot
    fig = px.scatter(candidates_df ,x='eps_ave', y='nat', color=np.log10(candidates_df['score']))
    return fig
    
def main():
    gen_plot()
    
if __name__ == "__main__":
    main()
    
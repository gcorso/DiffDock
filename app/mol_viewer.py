viewer_description = """
PDB viewer using 3Dmol.js  
If using please cite:  
3Dmol.js: Molecular visualization with WebGL, Nicholas Rego, David Koes , Bioinformatics, Volume 31, Issue 8, April 2015, Pages 1322–1324, https://doi.org/10.1093/bioinformatics/btu829

See also:
https://huggingface.co/blog/spaces_3dmoljs
https://huggingface.co/spaces/simonduerr/3dmol.js/blob/main/app.py
"""


def gen_3dmol_vis(pdb_text: str, sdf_text: str):
    x = (
            """<!DOCTYPE html>
            <html>
            <head>    
        <meta http-equiv="content-type" content="text/html; charset=UTF-8" />
        <style>
        body{
            font-family:sans-serif
        }
        .mol-container {
        width: 100%;
        height: 600px;
        position: relative;
        }
        .mol-container select{
            background-image:None;
        }
        </style>
         <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.3/jquery.min.js" integrity="sha512-STof4xm1wgkfm7heWqFJVn58Hm3EtS31XFaagaa8VMReCXAkQnJZ+jEy8PCC/iT18dFy95WcExNHFTqLyp72eQ==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        </head>
        <body>  
    
        <div id="container" class="mol-container"></div>
    
                <script>
                   let pdb = `""" + pdb_text + """`  
            
                   let sdf = `""" + sdf_text + """`

             $(document).ready(function () {
                let element = $("#container");
                let config = { backgroundColor: "white" };
                let viewer = $3Dmol.createViewer(element, config);
                
                viewer.addModel(pdb, "pdb", { format: "pdb" });
                pdb_model = viewer.getModel(0);
                // Cartoon view for protein
                // First argument is a selector
                pdb_model.setStyle({}, {cartoon: {color: "spectrum", opacity: 0.8}});
                
                viewer.addModel(sdf, "sdf", {format: "sdf"});
                sdf_model = viewer.getModel(1);
                // Stick view for SDF
                sdf_model.setStyle({}, {stick: {color: "red", radius: 0.2, opacity: 0.8}});
                
                viewer.zoomTo();
                viewer.render();
                // viewer.zoom(0.8, 2000);
              })
        </script>
        </body></html>"""
    )

    return f"""<iframe style="width: 100%; height: 600px" name="result" sandbox="allow-modals allow-forms 
    allow-scripts allow-same-origin allow-popups 
    allow-top-navigation-by-user-activation allow-downloads" allowfullscreen="" 
    allowpaymentrequest="" frameborder="0" srcdoc='{x}'></iframe>"""


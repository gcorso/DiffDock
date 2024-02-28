import collections
import datetime
import logging
import os
from typing import Tuple, Optional, Dict

import gradio as gr

import mol_viewer
import run_utils
from run_utils import PROJECT_URL, PROJECT_DIR, TEMP_DIR

DEFAULT_INFERENCE_ARGS = os.path.join(PROJECT_DIR, "default_inference_args.yaml")


def run_wrapper(protein_pdb_id, protein_file, ligand_smile, ligand_file, config_file, *args) -> Tuple[str, Optional[str], Optional[Dict], Optional[gr.Dropdown]]:

    if protein_pdb_id is not None and protein_file is None:
        protein_file_name = run_utils.download_pdb(protein_pdb_id, TEMP_DIR)
    else:
        protein_file_name = protein_file['name']

    if protein_file_name is None:
        return "Protein file is missing! Must provide a protein file in PDB format", None, None, None
    if ligand_file is None and ligand_smile is None:
        return "Ligand is missing! Must provide a ligand file in SDF format or SMILE string", None, None, None

    config_path = config_file['name'] if config_file else DEFAULT_INFERENCE_ARGS
    ligand_desc = ligand_file['name'] if ligand_file else ligand_smile
    output_file = run_utils.run_cli_command(
        protein_file_name, ligand_desc, config_path, *args,
    )

    message = f"Calculation completed at {datetime.datetime.now()}"

    view_selector_content = collections.OrderedDict()
    dropdown = None
    # print(f"Output file: {output_file}")
    if output_file:
        pdb_files, sdf_files = run_utils.process_zip_file(output_file)
        # print(f"PDB file: {pdb_files}")
        pdb_file = pdb_files[0] if pdb_files else None
        for sdf_file in sdf_files:
            confidence = sdf_file.get("confidence", None)
            # rank1 has no confidence
            if confidence is None:
                continue
            label = f"Rank {sdf_file['rank']}. Confidence {confidence:.2f}"
            pdb_text = pdb_file['content'] if pdb_file else None
            sdf_text = sdf_file['content']
            output_viz = "Output visualisation unavailable"
            if pdb_text:
                logging.debug(f"Creating 3D visualisation")
                output_viz = mol_viewer.gen_3dmol_vis(pdb_text, sdf_text)
            view_selector_content[label] = output_viz

        labels = list(view_selector_content.keys())
        init_value = labels[0] if labels else None
        dropdown = gr.Dropdown(interactive=True, label="Ranked samples",
                               choices=labels, value=init_value)

    return message, output_file, view_selector_content, dropdown


def update_view(view_selector_content, view_result_selector, default_str="Output visualisation unavailable"):
    if view_selector_content and view_result_selector:
        return view_selector_content.get(view_result_selector, default_str)
    return default_str


def run():

    with gr.Blocks(title="DiffDock Web") as demo:
        gr.Markdown("# DiffDock Web")
        gr.Markdown(f"""Run [DiffDock]({PROJECT_URL}) for a single protein and ligand.
        We have provided the most important inputs as UI elements.  """)
        with gr.Box():
            gr.Markdown("# Input")
            with gr.Row():
                with gr.Column():
                    gr.Markdown("## Protein")
                    protein_pdb_id = gr.Textbox(
                        placeholder="PDB Code or upload file below", label="Input PDB ID"
                    )
                    protein_pdb_file = gr.File(file_count="single", label="Input PDB File")
                with gr.Column():
                    gr.Markdown("## Ligand")
                    ligand_smile = gr.Textbox(
                        placeholder="Provide SMILES input or upload mol2/sdf file below",
                        label="SMILES string",
                    )
                    ligand_file = gr.File(file_count="single", label="Input Ligand", file_types=[".sdf", ".mol2"])

            with gr.Row():
                samples_per_complex = gr.Number(label="Samples Per Complex", value=10, minimum=1, maximum=100, precision=0, interactive=True)

            with gr.Row():
                with gr.Column():
                    config_instructions = f"""## Configuration (Optional)
                        Configuration file to be passed 
                        to [inference.py]({PROJECT_URL}/blob/main/inference.py). 
                        If this is provided, it must supply all necessary arguments.
                        If not provided, the [default configuration]({PROJECT_URL}/blob/main/app/default_inference_args.yml) will be used."""
                    gr.Markdown(config_instructions)

                    config_file = gr.File(label="Configuration (Optional, YML)", file_types=[".yml", ".yaml"], value=None,
                                          info="Additional arguments to pass to DiffDock.")

            with gr.Row():
                with gr.Column():
                    gr.Markdown("## Examples")
                    gr.Examples(
                        [
                            [
                                "6w70",
                                "examples/6w70.pdb",
                                "COc1ccc(cc1)n2c3c(c(n2)C(=O)N)CCN(C3=O)c4ccc(cc4)N5CCCCC5=O",
                                "examples/6w70_ligand.sdf",
                                10,
                                True
                            ],
                            [
                                "6moa",
                                "examples/6moa_protein_processed.pdb",
                                "",
                                "examples/6moa_ligand.sdf",
                                10,
                                True
                            ],
                            [
                                "",
                                "examples/6o5u_protein_processed.pdb",
                                "",
                                "examples/6o5u_ligand.sdf",
                                10,
                                True
                            ],
                            [
                                "",
                                "examples/6o5u_protein_processed.pdb",
                                "[NH3+]C[C@H]1O[C@H](O[C@@H]2[C@@H]([NH3+])C[C@H]([C@@H]([C@H]2O)O[C@H]2O[C@H](CO)[C@H]([C@@H]([C@H]2O)[NH3+])O)[NH3+])[C@@H]([C@H]([C@@H]1O)O)O",
                                "examples/6o5u_ligand.sdf",
                                10,
                                True
                            ],
                            [
                                "",
                                "examples/6o5u_protein_processed.pdb",
                                "",
                                "examples/6o5u_ligand.sdf",
                                10,
                                True
                            ],
                            [
                                "",
                                "examples/6ahs_protein_processed.pdb",
                                "",
                                "examples/6ahs_ligand.sdf",
                                10,
                                True
                            ],
                        ],
                        [protein_pdb_id, protein_pdb_file, ligand_smile, ligand_file, samples_per_complex],
                    )

        with gr.Row():
            run_btn = gr.Button("Run DiffDock")

        with gr.Box():
            gr.Markdown("# Output")
            with gr.Row():
                message = gr.Text(label="Run message", interactive=False)
            with gr.Row():
                output_file = gr.File(label="Output Files")
            with gr.Row():
                with gr.Column():
                    init_value = "DiffDock prediction visualization"
                    view_result_selector = gr.Dropdown(interactive=True, label="Ranked samples")
                    viewer = gr.HTML(value=init_value, label="Protein Viewer", show_label=True)

        with gr.Row():
            gr.Markdown("Many thanks to [Simon Duerr](https://huggingface.co/simonduerr), who created the "
                        "[original DiffDock web interface](https://huggingface.co/spaces/simonduerr/diffdock), "
                        "on which this interface is based.")

        view_selector_content = gr.Variable()

        _inputs = [protein_pdb_id, protein_pdb_file, ligand_smile, ligand_file, config_file]
        # See run_utils.py:ARG_ORDER for the order of these arguments
        _inputs += [samples_per_complex]
        _outputs = [message, output_file, view_selector_content, view_result_selector]
        run_btn.click(fn=run_wrapper, inputs=_inputs, outputs=_outputs, preprocess=False)

        view_result_selector.change(fn=update_view,
                                    inputs=[view_selector_content, view_result_selector],
                                    outputs=viewer)

    server_port = int(os.environ.get("GRADIO_SERVER_PORT", "7860"))
    demo.launch(server_name="0.0.0.0", server_port=server_port, share=False)


if __name__ == "__main__":
    run_utils.set_env_variables()
    run_utils.configure_logging()

    run()

import datetime
import io
import os
import shutil
import subprocess
import tempfile
import uuid

import logging
import zipfile
from typing import List, Dict

import requests

PROJECT_URL = "https://github.com/gcorso/DiffDock"

ARG_ORDER = ["samples_per_complex"]

APP_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(APP_DIR, ".."))
# Directory for semi-temporary files, ie PDB files downloaded from PDB
TEMP_DIR = os.path.join(APP_DIR, ".tmp")
os.makedirs(TEMP_DIR, exist_ok=True)


def set_env_variables():
    if "DiffDockDir" not in os.environ:
        work_dir = os.path.abspath(PROJECT_DIR)
        if os.path.exists(work_dir):
            os.environ["DiffDockDir"] = work_dir
        else:
            raise ValueError(f"DiffDockDir {work_dir} not found")

    if "LOG_LEVEL" not in os.environ:
        os.environ["LOG_LEVEL"] = "INFO"


def configure_logging(level=None):
    if level is None:
        level = getattr(logging, os.environ.get("LOG_LEVEL", "INFO"))

    # Note that this sets the universal logger,
    # which includes other libraries.
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] [%(filename)s:%(lineno)d] %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S %Z",
        handlers=[
            logging.StreamHandler(),  # Outputs logs to stderr by default
            # If you also want to log to a file, uncomment the following line:
            # logging.FileHandler('my_app.log', mode='a', encoding='utf-8')
        ],
    )


def kwargs_to_cli_args(**kwargs) -> List[str]:
    """
    Converts keyword arguments to a CLI argument string.
    Boolean kwargs are added as flags if True, and omitted if False.
    """
    cli_args = []
    for key, value in kwargs.items():
        if isinstance(value, bool):
            if value:
                cli_args.append(f"--{key}")
        else:
            if value is not None and str(value) != "":
                cli_args.append(f"--{key}={value}")

    return cli_args


def read_file_lines(fi_path: str, skip_remarks=True):
    with open(fi_path, "r") as fp:
        lines = fp.readlines()
        if skip_remarks:
            lines = list(filter(lambda x: not x.upper().startswith("REMARK"), lines))
    mol = "".join(lines)
    return mol


def run_cli_command(
    protein_path: str,
    ligand: str,
    config_path: str,
    *args,
    work_dir=None,
):
    if work_dir is None:
        work_dir = os.environ.get(
            "DiffDockDir", PROJECT_DIR
        )

    assert len(args) == len(ARG_ORDER), f'Expected {len(ARG_ORDER)} arguments, got {len(args)}'

    inference_log_level = os.environ.get("INFERENCE_LOG_LEVEL", os.environ.get("LOG_LEVEL", "WARNING"))

    all_arg_dict = {"protein_path": protein_path, "ligand": ligand, "config": config_path,
                    "no_final_step_noise": True, "loglevel": inference_log_level}
    for arg_name, arg_val in zip(ARG_ORDER, args):
        all_arg_dict[arg_name] = arg_val

    # Check device availability
    result = subprocess.run(
        ["python3", "utils/print_device.py"],
        cwd=work_dir,
        check=False,
        text=True,
        capture_output=True,
        env=os.environ,
    )
    logging.debug(f"Device check output:\n{result.stdout}")

    command = [
        "python3",
        "inference.py"]

    command += kwargs_to_cli_args(**all_arg_dict)

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = temp_dir
        command.append(f"--out_dir={temp_dir_path}")

        # Convert command list to string for printing
        command_str = " ".join(command)
        logging.info(f"Executing command: {command_str}")

        # Running the command
        try:
            # Option to skip running for UI-dev purposes
            skip_running = os.environ.get("__SKIP_RUNNING", "false").lower() == "true"
            if not skip_running:
                result = subprocess.run(
                    command,
                    cwd=work_dir,
                    check=False,
                    text=True,
                    capture_output=True,
                )
                logging.debug(f"Command output:\n{result.stdout}")
                full_output = f"Standard out:\n{result.stdout}"
                if result.stderr:
                    # Skip progress bar lines
                    stderr_lines = result.stderr.split("\n")
                    stderr_lines = filter(lambda x: "%|" not in x, stderr_lines)
                    stderr_text = "\n".join(stderr_lines)
                    logging.error(f"Command error:\n{stderr_text}")
                    full_output += f"\nStandard error:\n{stderr_text}"

                with open(f"{temp_dir_path}/output.log", "w") as log_file:
                    log_file.write(full_output)

            else:
                logging.debug("Skipping command execution")
                artificial_output_dir = os.path.join(TEMP_DIR, "artificial_output")
                os.makedirs(artificial_output_dir, exist_ok=True)
                shutil.copy(protein_path, os.path.join(artificial_output_dir, "protein.pdb"))
                shutil.copy(ligand, os.path.join(artificial_output_dir, "rank1.sdf"))
                shutil.copy(ligand, os.path.join(artificial_output_dir, "rank1_confidence-0.10.sdf"))

        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred while executing the command: {e}")

        # Copy the input protein into the output directory
        sub_dirs = [os.path.join(temp_dir_path, x) for x in os.listdir(temp_dir_path)]
        sub_dirs = list(filter(lambda x: os.path.isdir(x), sub_dirs))
        logging.debug(f"Output Subdirectories: {sub_dirs}")
        if len(sub_dirs) == 1:
            sub_dir = sub_dirs[0]
            # Copy the input protein from the input to the output
            trg_protein_path = os.path.join(sub_dir, os.path.basename(protein_path))
            shutil.copy(protein_path, trg_protein_path)

        # Zip the output directory
        # Generate a unique filename using a timestamp and a UUID
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        uuid_tag = str(uuid.uuid4())[0:8]
        unique_filename = f"diffdock_output_{timestamp}_{uuid_tag}"
        zip_base_name = os.path.join("tmp", unique_filename)

        logging.debug(f"About to zip directory '{temp_dir}' to {unique_filename}")

        full_zip_path = shutil.make_archive(zip_base_name, "zip", temp_dir)

        logging.debug(f"Directory '{temp_dir}' zipped to {unique_filename}'")

    return full_zip_path


def parse_ligand_filename(filename: str) -> Dict:
    """
    Parses an sdf filename to extract information.
    """
    if not filename.endswith(".sdf"):
        return {}

    basename = os.path.basename(filename).replace(".sdf", "")
    tokens = basename.split("_")
    rank = tokens[0]
    rank = int(rank.replace("rank", ""))
    if len(tokens) == 1:
        return {"filename": basename, "rank": rank, "confidence": None}

    con_str = tokens[1]
    conf_val = float(con_str.replace("confidence", ""))

    return {"filename": basename, "rank": rank, "confidence": conf_val}


def process_zip_file(zip_path: str):
    pdb_file = []
    sdf_files = []
    with zipfile.ZipFile(open(zip_path, "rb")) as my_zip_file:
        for filename in my_zip_file.namelist():
            # print(f"Processing file {filename}")
            if filename.endswith("/"):
                continue

            if filename.endswith(".pdb"):
                content = my_zip_file.read(filename).decode("utf-8")
                pdb_file.append({"path": filename, "content": content})

            if filename.endswith(".sdf"):
                info = parse_ligand_filename(filename)
                info["content"] = my_zip_file.read(filename).decode("utf-8")
                info["path"] = filename
                sdf_files.append(info)

    sdf_files = sorted(sdf_files, key=lambda x: x.get("rank", 1_000))

    return pdb_file, sdf_files


def download_pdb(pdb_code: str, work_dir: str):
    pdb_url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    pdb_path = os.path.join(work_dir, f"{pdb_code}.pdb")
    if not os.path.exists(pdb_path):
        logging.debug(f"Downloading PDB file for {pdb_code} from {pdb_url}")
        response = requests.get(pdb_url, allow_redirects=True)
        if response.status_code == 200:
            with open(pdb_path, "w") as pdb_file:
                pdb_file.write(response.text)
        else:
            logging.error(f"Failed to download PDB file for {pdb_code} from {pdb_url}")
            pdb_path = None

    else:
        logging.info(f"PDB file for {pdb_code} already exists at {pdb_path}")

    return pdb_path


def test_run_cli():
    # Testing code
    set_env_variables()
    configure_logging()

    work_dir = os.path.abspath(PROJECT_DIR)
    os.environ["DiffDockDir"] = work_dir
    protein_path = os.path.join(work_dir, "data", "3dpf", "3dpf_protein.pdb")
    ligand = os.path.join(work_dir, "data", "3dpf", "3dpf_ligand.sdf")
    config_file = os.path.join(APP_DIR, "default_inference_args.yaml")

    run_cli_command(
        protein_path,
        ligand,
        config_file,
        10,
        False,
        True,
        None
    )

import json
from argparse import ArgumentParser
from diffdock_protocol import DiffDockProtocol
from diffdock_api import DiffDockApi

# CLI version of DiffDockApi. Implemented for testing

def main(args):
    print(f"Running DiffDock CLI with {args}")
    with open(args.input, 'r', encoding='utf8') as inf:
        inputJson = inf.read()
    request = DiffDockProtocol.Request.from_json(inputJson)
    response = DiffDockApi.run_diffdock(request)
    with open(args.output, 'w', encoding='utf8') as outf:
        json.dump(response, outf, default=lambda x: x.__dict__)
    print(f"Done with DiffDock CLI")

if __name__ == "__main__":
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--input', type=str, help='json file with diffdock request', required=True)
    arg_parser.add_argument('--output', type=str, help='json file with diffdock response', required=True)
    args = arg_parser.parse_args()    
    main(args)

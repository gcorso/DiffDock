import json
from dataclasses import asdict, astuple, dataclass, field
from typing import List, Any
from collections import namedtuple

class DiffDockProtocol:
    MessageType = namedtuple('MessageType', ['INVALID', 'ERROR', 'RESULTS', 'STATUS'])('invalid', 'error', 'results', 'status')

    @dataclass
    class Protein:
        name: str
        data: str

        def __iter__(self):
            return iter(astuple(self))

    @dataclass
    class Ligand:
        name: str
        data: str

        def __iter__(self):
            return iter(astuple(self))

    @dataclass
    class Request:
        proteins: List['DiffDockProtocol.Protein'] = field(default_factory=list)
        ligands: List['DiffDockProtocol.Ligand'] = field(default_factory=list)
        samples_per_complex: int = 10
        add_hs: bool = True
        keep_hs: bool = False
        keep_src_3d: bool = False

        @staticmethod
        def from_json(data: str):
            return DiffDockProtocol.Request.from_dict(json.loads(data))

        @staticmethod
        def from_dict(data: dict):
            req = DiffDockProtocol.Request(
                proteins=[DiffDockProtocol.Protein(**p) for p in data['proteins']],
                ligands=[DiffDockProtocol.Ligand(**l) for l in data['ligands']],
            )
            # Need to preserve defaults for these if not in request
            if data.get('samples_per_complex') is not None: req.samples_per_complex = data['samples_per_complex']
            if data.get('add_hs') is not None: req.add_hs = data['add_hs']
            if data.get('keep_hs') is not None: req.keep_hs = data['keep_hs']
            if data.get('keep_src_3d') is not None: req.keep_src_3d = data['keep_src_3d']
            return req
            
    @dataclass
    class Pose:
        filename: str
        sdf: str
        rank: int
        confidence: float

        def __iter__(self):
            return iter(astuple(self))

    @dataclass
    class Result:
        proteinName: str
        ligandName: str
        poses: List['DiffDockProtocol.Pose'] = field(default_factory=list)
        error: str = ''

        @staticmethod
        def from_dict(data: dict):
            return DiffDockProtocol.Result(
                proteinName=data['proteinName'],
                ligandName=data['ligandName'],
                poses=[DiffDockProtocol.Pose(**p) for p in data['poses']],
                error=data['error']
            )

    @dataclass
    class Response:
        messageType: 'DiffDockProtocol.MessageType' = 'invalid'
        results: List['DiffDockProtocol.Result'] = field(default_factory=list)
        message: str = ''

        @staticmethod
        def makeError(message: str=""):
            return DiffDockProtocol.Response(messageType=DiffDockProtocol.MessageType.ERROR, message=message)

        @staticmethod
        def makeStatus(message: str=""):
            return DiffDockProtocol.Response(messageType=DiffDockProtocol.MessageType.STATUS, message=message)

        @staticmethod
        def makeResults(results: List['DiffDockProtocol.Result']=[]):
            return DiffDockProtocol.Response(messageType=DiffDockProtocol.MessageType.RESULTS, results=results)

        @staticmethod
        def from_json(data: str):
            return DiffDockProtocol.Response.from_dict(json.loads(data))

        @staticmethod
        def from_dict(data: dict):
            messageType = data['messageType']
            if messageType not in DiffDockProtocol.MessageType:
                messageType = DiffDockProtocol.MessageType.INVALID
            return DiffDockProtocol.Response(
                messageType=messageType,
                results=[DiffDockProtocol.Result.from_dict(r) for r in data['results']],
                message=data['message']
            )

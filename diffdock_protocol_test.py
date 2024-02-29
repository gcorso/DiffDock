import json

from diffdock_protocol import DiffDockProtocol


def main():
    print("Testing diffdock_protocol")
    print("Running request tests...")
    test_requests()
    print("Running response tests...")
    test_responses()
    print("All tests passed!")

def test_requests():
    protein1 = DiffDockProtocol.Protein("protein1", "data1")
    protein2 = DiffDockProtocol.Protein("protein2", "data2")
    protein3 = DiffDockProtocol.Protein("protein3", "data3")
    protein4 = DiffDockProtocol.Protein("protein4", "data4")
    ligand1 = DiffDockProtocol.Ligand("ligand1", "data1")
    ligand2 = DiffDockProtocol.Ligand("ligand2", "data2")
    ligand3 = DiffDockProtocol.Ligand("ligand3", "data3")
    ligand4 = DiffDockProtocol.Ligand("ligand4", "data4")

    request = DiffDockProtocol.Request()
    check_request(request)
    request.proteins.append(DiffDockProtocol.Protein("protein1", "data1"))
    request.ligands.append(DiffDockProtocol.Ligand("ligand1", "data1"))
    request.add_hs = True
    request.keep_hs = False
    request.keep_src_3d = False
    request.samples_per_complex = 10
    check_request(request)

    request = DiffDockProtocol.Request()
    request.proteins.append(DiffDockProtocol.Protein("protein1", "data1"))
    request.proteins.append(DiffDockProtocol.Protein("protein2", "data2"))
    request.ligands.append(DiffDockProtocol.Ligand("ligand1", "data1"))
    request.ligands.append(DiffDockProtocol.Ligand("ligand2", "data2"))
    check_request(request)

    request = DiffDockProtocol.Request([protein1, protein2, protein3], [ligand1, ligand2, ligand3], 5, False, True, True)
    check_request(request)

    request = DiffDockProtocol.Request([protein2, protein3], [ligand1, ligand3], 6)
    check_request(request)
    request = DiffDockProtocol.Request([protein1, protein3], [ligand2], 7, False)
    check_request(request)

def test_responses():
    response = DiffDockProtocol.Response()
    check_response(response)
    response.messageType = DiffDockProtocol.MessageType.RESULTS
    response.results.append(DiffDockProtocol.Result("protein1", "ligand1", [DiffDockProtocol.Pose("filename1", "sdf1", 1, 0.9)]))
    response.error = ""
    check_response(response)

    response = DiffDockProtocol.Response.makeError("There was a problem")
    check_response(response)

    response = DiffDockProtocol.Response.makeError("There was a problem2")
    check_response(response)
    
    poses = [
        DiffDockProtocol.Pose("file1", "data1", 1, 1.0),
        DiffDockProtocol.Pose("file2", "data2", 2, 0.0),
        DiffDockProtocol.Pose("file3", "data3", 3, -1.0),
    ]
    result = DiffDockProtocol.Result("protein1", "ligand1", poses)
    results = [result]
    response = DiffDockProtocol.Response.makeResults(results)
    check_response(response)

    result = DiffDockProtocol.Result("protein2", "ligand2", error="This complex failed")
    response = DiffDockProtocol.Response.makeResults([result])
    check_response(response)

    response = DiffDockProtocol.Response(messageType="fail")
    response2 = DiffDockProtocol.Response(messageType="invalid")
    check_response(response, AssertionError)
    check_responses(serializeAndDeserializeResponse(response), response2)

# Request helpers
def check_requests(request1, request2):
    assert len(request1.proteins) == len(request2.proteins)
    assert len(request1.ligands) == len(request2.ligands)
    assert request1.add_hs == request2.add_hs
    assert request1.keep_hs == request2.keep_hs
    assert request1.keep_src_3d == request2.keep_src_3d
    assert request1.samples_per_complex == request2.samples_per_complex

    for i in range(len(request1.proteins)):
        protein1name = request1.proteins[i].name
        protein2name = request2.proteins[i].name
        assert protein1name == protein2name
        assert request1.proteins[i].data == request2.proteins[i].data

    for i in range(len(request1.ligands)):
        assert request1.ligands[i].name == request2.ligands[i].name
        assert request1.ligands[i].data == request2.ligands[i].data

def check_request(request):
    request2 = serializeAndDeserializeRequest(request)
    check_requests(request, request2)

    print(f"Request test passed")

def serializeAndDeserializeRequest(request):
    j = json.dumps(request, default=lambda x: x.__dict__)
    request2 = DiffDockProtocol.Request.from_json(j)
    return request2

# Response helpers
def check_responses(response1, response2):
    assert response1.messageType == response2.messageType
    assert len(response1.results) == len(response2.results)
    assert response1.message == response2.message

    for i in range(len(response1.results)):
        result1 = response1.results[i]
        result2 = response2.results[i]

        assert result1.proteinName == result2.proteinName
        assert result1.ligandName == result2.ligandName
        assert result1.error == result2.error
        assert len(result1.poses) == len(result2.poses)
        for j in range(len(result1.poses)):
            pose1 = result1.poses[j]
            pose2 = result2.poses[j]

            assert pose1.filename == pose2.filename
            assert pose1.sdf == pose2.sdf
            assert pose1.rank == pose2.rank
            assert pose1.confidence == pose2.confidence

def check_response(response, exceptionExpected=None):
    response2 = serializeAndDeserializeResponse(response)
    haveException = False
    try:
        check_responses(response, response2)
    except Exception as ex:
        if not exceptionExpected or not isinstance(ex, exceptionExpected):
            raise
        haveException = True
    if exceptionExpected is not None:
        assert haveException

    print(f"Response test passed")


def serializeAndDeserializeResponse(response):
    j = json.dumps(response, default=lambda x: x.__dict__)
    response2 = DiffDockProtocol.Response.from_json(j)
    return response2

# main
if __name__ == "__main__":
    main()

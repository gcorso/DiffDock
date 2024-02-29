import asyncio
import json
import logging
import re
from websockets.server import serve
import concurrent.futures
from argparse import ArgumentParser

from diffdock_protocol import DiffDockProtocol
from diffdock_api import DiffDockApi

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('diffdock_service')

#
# This is a websocket server for servicing DiffDock requests.
# This listens on a port and places incoming json requests into a queue, serviced by a number of workers.
#

async def handleRequest(websocket, queue):
    """
    Receive a request and put it into the queue.
    Make sure we keep the websocket so we can send the response when docking finishes.
    """
    message = await websocket.recv()
    log.info("Received a request.")
    requestId, cmdName, requestData = extractWsAppRequest(message)
    requestContext = websocket, requestId

    try:
        requestObj = DiffDockProtocol.Request.from_json(requestData)
    except:
        log.exception(f"Rejected a request")
        await sendError(requestContext, f"Invalid request")
        return

    log.info("handleRequest running a request...")
    try:
        # Put the request into the queue with a condition to resolve when it finishes.
        # Otherwise the websocket will close early
        condition = asyncio.Condition()
        await queue.put((requestObj, requestContext, condition))
        await sendStatus(requestContext, "Request received")
        async with condition:
            await condition.wait()
    except Exception as ex:
        log.exception(f"handleRequest encountered exception")
        await sendError(requestContext, f"Service exception: {ex}")


async def queueWorker(queue):
    """
    Handle docking requests from the queue
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        while True:
            # Wait for an item from the queue
            requestObj, requestContext, condition = await queue.get()

            log.info("worker running a request...")
            try:
                # Since this is async, but diffdock is not, it seems it is tricky.
                # run_in_executor seemed to do the trick.
                # We want to make sure we don't block other connections from putting their requests
                # into the queue.
                await sendStatus(requestContext, "Working request...")
                loop = asyncio.get_event_loop()
                response = await loop.run_in_executor(executor, DiffDockApi.run_diffdock, requestObj)
                log.info("queueWorker resolved response")
                await sendResults(requestContext, response)
                log.info("queueWorker sent response.")
            except Exception as ex:
                logger.exception(f"queueWorker encountered exception")
                await sendError(requestContext, f"Service exception: {ex}")

            # Notify the queue that the item has been processed
            queue.task_done()

            # Resolve the condition for the request handler
            async with condition:
                condition.notify_all()


async def main(host="localhost", port=9002, max_size=2**24, worker_count=5):
    queue = asyncio.Queue()
    workers = [asyncio.create_task(queueWorker(queue)) for _ in range(worker_count)]

    async def handlerWrapper(websocket, path):
        await handleRequest(websocket, queue)

    log.info(f"DiffDock service starting websocket listener at {host}:{port}...")
    start_server = serve(handlerWrapper, host, port, max_size=max_size)

    log.info(f"DiffDock ready.")
    try:
        await start_server
        await asyncio.Future()
    except KeyboardInterrupt:
        log.info("Received KeyboardInterrupt")

    # SHUTDOWN
    try:
        log.info("DiffDock service shutting down...let this finish naturally if you can")
        server.close()
        loop.run_until_complete(server.wait_closed())
        for worker in workers:
            worker.cancel()
    except:
        log.info("DiffDock service caught shutting down, will force termination")
    finally:
        log.info("DiffDock service finished.")
        os._exit(0)
        

arg_parser = ArgumentParser()
arg_parser.add_argument('--port', type=int, default=9002, help='Port to listen on')
args = arg_parser.parse_args()
log.info("Starting inference service...")
asyncio.run(main(port=args.port))

## Helpers
def extractWsAppRequest(req):
    """
    Bioleap WsApps requests look like `<requestid> <cmdname> <request data>`
    """
    match = re.match("(#[a-z-_]+\d+) ([a-z-_]+) (\{.*\})", req)
    if match:
        requestId, cmdName, data = match.groups()
        return requestId, cmdName, data
    else:
        return None, None, req


async def sendPacket(requestContext, responseObj):
    responseData = json.dumps(responseObj, default=lambda x: x.__dict__)
    websocket, requestId = requestContext

    if not requestId:
        return await websocket.send(responseData)

    # We have a Bioleap WsApp request id
    responseCmd = None
    messageType = responseObj.messageType
    complete = False
    if messageType == DiffDockProtocol.MessageType.ERROR:
        responseCmd = "diffdock-error"
        complete = True
    elif messageType == DiffDockProtocol.MessageType.RESULTS:
        responseCmd = "diffdock-results"
        complete = True
    elif messageType == DiffDockProtocol.MessageType.STATUS:
        responseCmd = "diffdock-status"

    payloadParts = [x for x in [requestId, responseCmd, responseData] if x]

    payload = ' '.join(payloadParts)
    log.info(f'Sending {payload}')
    await websocket.send(payload)
    if complete:
        # For Bioleap WsApps we need to end with a `complete` message
        log.info(f'Sending complete')
        completePayload = ' '.join([requestId, 'completed'])
        return await websocket.send(completePayload)

async def sendStatus(requestContext, status):
    responseObj = DiffDockProtocol.Response.makeStatus(status)
    return await sendPacket(requestContext, responseObj)

async def sendError(requestContext, error):
    responseObj = DiffDockProtocol.Response.makeError(error)
    return await sendPacket(requestContext, responseObj)

async def sendResults(requestContext, responseObj):
    return await sendPacket(requestContext, responseObj)

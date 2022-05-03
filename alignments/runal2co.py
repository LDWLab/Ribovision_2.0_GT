

def executeAl2co(alnstring):
    from subprocess import Popen, PIPE
    from os import remove, path
    import datetime, re

    now = datetime.datetime.now()
    fileNameSuffix = "_" + str(now.year) + "_" + str(now.month) + "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) + "_" + str(now.second) + "_" + str(now.microsecond)
    clustalName = f"./static/aln2co{fileNameSuffix}.aln"

    if path.isfile(clustalName):
        remove(clustalName)

    fh = open(clustalName, "w")
    fh.write(alnstring.replace('CLUSTAL X', 'CLUSTAL W'))
    fh.close()

    pipe = Popen(f"al2co -i {clustalName} -c 2", stdout=PIPE, shell=True)
    output = pipe.communicate()[0]

    if len(output.decode("ascii")) <= 0:
        remove(clustalName)
        return None

    al2coOutput = output.decode("ascii")
    al2coResult = dict()
    for posResult in al2coOutput.split('\n'):
        if re.match(r'^\*', posResult):
            break
        al2coResult[int(posResult.split()[0])] = float(posResult.split()[2])

    remove(clustalName)
    return al2coResult
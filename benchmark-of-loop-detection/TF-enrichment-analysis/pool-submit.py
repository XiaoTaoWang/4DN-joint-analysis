import subprocess, glob, os

queue = glob.glob('bigwigs/*bw')
loop_fil = 'H1ESC.union-loops.with-cluster-info.bedpe'
for bw in queue:
    TF = os.path.split(bw)[1].split('.')[1]

    Indicator = '{0}.completed'.format(TF)
    lockFile = '{0}.lock'.format(TF)

    if os.path.exists(Indicator):
        continue
    
    if os.path.exists(lockFile):
        continue
    
    lock = open(lockFile, 'w')
    lock.close()

    outfil = '{0}-{1}.matrices-cache.pkl'.format('H1', TF)
    command = ['python', 'extract-matrices.py', loop_fil, bw, outfil]
    subprocess.check_call(' '.join(command), shell=True)

    completed = open(Indicator, 'wb')
    completed.close()

    os.remove(lockFile)
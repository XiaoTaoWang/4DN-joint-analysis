import os, subprocess, cooler

queue = [
    ('H1ESC', 'Hi-C', 'H1ESC_Hi-C_FA_DSG_DpnII', 'HIC-H1ESC-DpnII-allReps.4DNFI82R42AD.mcool'),
    ('H1ESC', 'Micro-C', 'H1ESC_MicroC', 'MicroC-H1ESC-MNase-allReps.mcool'),
    ('H1ESC', 'PLAC-Seq', 'H1ESC_PLAC-Seq', 'PLAC-H1ESC-H3K4me3-allReps.mcool'),
    ('H1ESC', 'ChIA-PET_CTCF', 'H1ESC_ChIA-PET_CTCF', 'H1ESC.CTCF-ChIAPET.4DNFINMHXGVQ.mcool'),
    ('H1ESC', 'ChIA-PET_Pol2', 'H1ESC_ChIA-PET_Pol2', 'H1ESC.RNAPII-ChIAPET.4DNFIO8YJ5JA.mcool'),
    ('HFFc6', 'Hi-C', 'HFFc6_Hi-C_FA_DSG_DpnII', 'HIC-HFFc6-DpnII-allReps.4DNFIAVXXO55.mcool'),
    ('HFFc6', 'Micro-C', 'HFFc6_MicroC', 'MicroC-HFFc6-MNase-allReps.mcool'),
    ('HFFc6', 'PLAC-Seq', 'HFFc6_PLAC-Seq', 'PLAC-HFFc6-H3K4me3-allReps.mcool'),
    ('HFFc6', 'ChIA-PET_CTCF', 'HFFc6_ChIA-PET_CTCF', 'HFFc6.CTCF-ChIAPET.4DNFIG2ILS39.mcool'),
    ('HFFc6', 'ChIA-PET_Pol2', 'HFFc6_ChIA-PET_Pol2', 'HFFc6.RNAPII-ChIAPET.4DNFIIASUSSX.mcool')
]

resolutions = [2000]
for res in resolutions:
    for cell, platform, label, cool_fil in queue:
        uri = '{0}::resolutions/{1}'.format(cool_fil, res)

        Indicator = '{0}.{1}.completed'.format(label, res)
        lockFile = '{0}.{1}.lock'.format(label, res)
        if os.path.exists(Indicator):
            continue
        
        if os.path.exists(lockFile):
            continue
        
        lock = open(lockFile, 'w')
        lock.close()

        outfil = '{0}.Enhancer_backgrounds.{1}kb.pkl'.format(label, res//1000)
        command = ['python', 'extract-loops-highres.py', uri, cell, outfil]
        subprocess.check_call(' '.join(command), shell=True)
        
        completed = open(Indicator, 'wb')
        completed.close()

        os.remove(lockFile)
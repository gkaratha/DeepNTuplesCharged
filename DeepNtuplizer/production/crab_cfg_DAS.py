from CRABClient.UserUtilities import config, ClientException
from CRABAPI.RawCommand import crabCommand
import datetime
import os.path, subprocess


def list_of_files(path):
  files= subprocess.check_output(["ls", "/eos/cms/"+path]).splitlines()
  outfiles=''
  iline =0 
  for line in files:
    if b"root" in line:
       if iline<len(files):
          outfiles+="'"+path+(line.decode())+"',"
       else:
          outfiles+="'"+path+(line.decode())+"'"
    iline+=1
  return outfiles

def submit(config):
    crabCommand('submit', config = config)


templ_sub="from WMCore.Configuration import Configuration\n"
templ_sub += "config = Configuration()\n"
templ_sub += "config.section_('General')\n"
templ_sub +="config.General.transferOutputs = True\n"
###### name local
templ_sub +="config.General.workArea = 'crab_projects/PFC_TTbar'\n"
templ_sub +="config.section_('JobType')\n"
templ_sub +="config.JobType.pluginName = 'Analysis'\n"
####### Here put the code
templ_sub +="config.JobType.psetName = 'DeepNtuplizer_pfc2.py'\n"
######
templ_sub +="config.JobType.allowUndistributedCMSSW = True\n"
templ_sub +="config.JobType.maxMemoryMB = 3500\n"
templ_sub +='config.JobType.inputFiles = ["../python/QGL_cmssw8020_v2.db"]\n'
templ_sub +="config.section_('Data')\n"
templ_sub +="config.Data.splitting = 'FileBased'\n"
templ_sub +="config.Data.unitsPerJob = 1\n"
templ_sub +="config.Data.inputDBS = 'global'\n"
templ_sub +="config.Data.publication = False\n"
###### Name in eos
templ_sub +="config.Data.outputDatasetTag = 'PFC_TTbar'\n"
templ_sub +="config.section_('Site')\n"
templ_sub +="config.Site.storageSite = 'T2_CH_CERN'\n"
##### output folder
templ_sub +="config.Data.outLFNDirBase = '/store/group/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v2/'\n"
##### input
templ_sub += "config.Data.inputDataset = '/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM'"
### DY2e4jets
# /DYto2E-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM
### QCD
# /QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM
# /QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26_ext1-v2/MINIAODSIM
### Wjets
# /WtoENu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM
### TTbar
# /TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM

with open("to_sub.py","w") as txt:
   txt.write(templ_sub)
txt.close()
os.system("crab submit -c to_sub.py")


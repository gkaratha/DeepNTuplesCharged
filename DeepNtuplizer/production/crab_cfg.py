from CRABClient.UserUtilities import config, ClientException
from CRABAPI.RawCommand import crabCommand
import datetime
import os.path, subprocess


def list_of_files(path):
  files= subprocess.check_output(["ls", "/eos/cms/"+path]).splitlines()
  outfiles=''
  for line in files:
    if b"root" in line:
       outfiles+=path+(line.decode())+"','"
  return outfiles

def submit(config):
    crabCommand('submit', config = config)


templ_sub="from WMCore.Configuration import Configuration\n"
templ_sub += "config = Configuration()\n"
templ_sub += "config.section_('General')\n"
templ_sub +="config.General.transferOutputs = True\n"
templ_sub +="config.General.workArea = 'crab_projects/ChargedJets_vsPuppi_v2_25_01_25'\n"
templ_sub +="config.section_('JobType')\n"
templ_sub +="config.JobType.pluginName = 'Analysis'\n"
templ_sub +="config.JobType.psetName = 'runGenBanalysis.py'\n"
templ_sub +="config.JobType.allowUndistributedCMSSW = True\n"
templ_sub +="config.JobType.maxMemoryMB = 3500\n"
templ_sub +="config.section_('Data')\n"
templ_sub +="config.Data.splitting = 'FileBased'\n"
templ_sub +="config.Data.unitsPerJob = 10\n"
templ_sub +="config.Data.inputDBS = 'global'\n"
templ_sub +="config.Data.publication = False\n"
templ_sub +="config.Data.outputDatasetTag = 'ChargedJets_vsPuppi_v2_25_01_25'\n"
templ_sub +="config.section_('Site')\n"
templ_sub +="config.Site.storageSite = 'T2_CH_CERN'\n"
templ_sub +="config.Data.outLFNDirBase = '/store/group/cmst3/user/gkaratha'\n"
templ_sub +="config.Data.userInputFiles = ["+list_of_files('/store/cmst3/group/softJets/gkaratha/chain_m70_dm20_cfgRun24_133X_Run2024_test_10172024/Mini/')+"]\n"
with open("to_sub.py","w") as txt:
   txt.write(templ_sub)
txt.close()
os.system("crab submit -c to_sub.py")


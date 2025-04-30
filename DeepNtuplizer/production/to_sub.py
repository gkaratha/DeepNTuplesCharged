from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects/PFC_TTbar'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'DeepNtuplizer_pfc2.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 3500
config.JobType.inputFiles = ["../python/QGL_cmssw8020_v2.db"]
config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'PFC_TTbar'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/group/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v2/'
config.Data.inputDataset = '/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM'
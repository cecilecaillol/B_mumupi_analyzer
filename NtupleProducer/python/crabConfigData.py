from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Bmmp_BPH5Run2018D'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_ntuplizer_data_cfg.py'
config.JobType.maxJobRuntimeMin = 400
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/ParkingBPH5/Run2018D-05May2019promptD-v1/MINIAOD'
config.Data.allowNonValidInputDataset = True

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100 # Run2018D 110010
config.Data.totalUnits = -1
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
#config.Data.outLFNDirBase = '/store/user/wangjian/'
config.Data.outLFNDirBase = '/store/user/ccaillol/'

config.Site.storageSite = 'T2_CH_CERN'

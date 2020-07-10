from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Bmmp_signal'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_ntuplizer_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 400

config.Data.inputDataset = '/BptoMMP_bbbar/wangjian-RunIIAutumn18MiniAOD-102X-c64c34953e011388ee7948c74fdb7fa0/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/ccaillol/'
#config.Data.outLFNDirBase = '/store/user/wangjian/'

config.Site.storageSite = 'T2_CH_CERN'


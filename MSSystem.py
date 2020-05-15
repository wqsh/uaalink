# -*- coding: utf-8 -*-

VALUE_MAX_SCAN = 500000
VALUE_APP_SCAN = 100000


EXPIRATION_TIME = {'Year': 2020, 'Month': 12, 'Day': 31}


INFO_TO_USER_Staff = (
    '\n[UUALink] UUALink is expired! Please send e-mail to cliu126@126.com for the new version.',
    '\n[UUALink] Warning! The current license will expired in 7 days. Please send e-mail to cliu126@126.com for the new version.',
    '\n[UUALink] Starting...',
    '\n[UUALink] Finished!',
    '\n[UUALink] Writing config file...',)

INFO_TO_USER_TaskRead = (
    '\n[UUALink] Reading MS1 file...',
    '\n[UUALink] Reading MS2 file...',
    '\n[UUALink] Reading identification results: ',)


INFO_TO_USER_TaskDraw = (
    '\n[UUALink] Drawing figures...',)

INFO_TO_USER_TaskExport = (
    '\n[UUALink] Exporting...',)

INFO_TO_USER_TaskXtract = (

    '\n[UUALink] Xtracting...',
    '\n[UUALink] ms1 and ms2 files are existed...',)


INFO_TO_USER_TaskProfile = (

    '\n[UUALink] Creating index file for ms1...',
    '\n[UUALink] Reconstructing chromatograms...',
    '\n[UUALink] #Chromatograms: ',)


FILENAME_EXPORT = (
    'INFO_MS1.txt',
    'INFO_MS2.txt',
    'INFO_Cycle.txt',
    'INFO_ID.txt',
    'INFO_Mass_Deviation.txt',
    'INFO_Chromatography.txt',)

FILENAME_DRAW = (
    'Fig_Histogram_IonInjectionTime_MS1.pdf',
    'Fig_Histogram_IonInjectionTime_MS2.pdf',
    'Fig_Histogram_ElutingTime.pdf',)


TYPE_IDENTIFICATION_RESULT = {'UUALink': 0, }

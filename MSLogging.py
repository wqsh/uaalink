import logging
import sys


def logGetError(info):
    print(info)

    logging.basicConfig(filename='UltraPlot.log',
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.error(info)
    sys.exit(0)
import os
import sys
import copy
import glob
import random
import contextlib

try:
    from urllib.parse import urlencode

except ImportError:
    from urllib import urlencode

try:
    from urllib.request import urlopen

except ImportError:
    from urllib2 import urlopen


def make_tiny(url):
  request_url = ('http://tinyurl.com/api-create.php?' +
  urlencode({'url':url}))

  with contextlib.closing(urlopen(request_url)) as response:
    return response.read().decode('utf-8').split('//')[1].replace('tinyurl.com/', '')


if __name__ == '__main__':
  turl = make_tiny('http://viewer.legacysurvey.org/?ra=334.7589&dec=0.7745&zoom=16&layer=dr8')

  print(turl)

  print('\n\nDone.\n\n')


#import xml.etree.cElementTree as etree
import xml.etree.ElementTree as ET
import itertools
from pprint import pprint

xmlstr = '/media/ngs/data/ncbi_expression/test.xml'
xmlstr = '/media/ngs/data/ncbi_expression/PRJEB4337.108.expression.xml'

with open(xmlstr) as f:
    it = itertools.chain('<root>', f, '</root>')
    root = ET.fromstringlist(it)

# Do something with `root`
root.find('.//tag3')



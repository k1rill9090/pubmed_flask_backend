from datetime import datetime
import sys
import os
import re
# print(os.getcwd())
# dt = datetime.strptime('2023-05-13', '%Y-%m-%d').strftime('%Y/%m/%d')
# print(type(dt))
# a = '     cov'
# b = a.strip()
# c = a
# print(b)
# print(c)

ss = 'fdfd  gf hg '
# match = re.findall(r'"([a-zA-Z\s\d-]+?)"', ss)
match = re.search(r'[^a-zA-Z\s\d-]+?', ss)
print(match)

# def doit(text):
#   import re
#   matches = re.findall(r'"([a-zA-Z\s\d-]+?)"',text)
# #   matches is now ['String 1', 'String 2', 'String3']
#   return ",".join(matches)
# # a = doit('Regex should return "String 1" or "String 2" or "String3" ')
# a = doit('"fdfd" "2fdg" "fd" "11"')
# print(a)
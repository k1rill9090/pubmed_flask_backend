# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 20:31:29 2023

@author: PC
"""

def FormStr(S,BegD,EndD):
    Words = S.split()
#    print(Words)
#    print(len(Words))
    TERM1 = ''
    if   len(Words) == 1:
    # TERM = '((((calcification[All Fields]) AND breast[All Fields]) AND cancer[All Fields]) AND pattern[ABSTRACT]) AND ("1970/01/01"[PubDate] : "2023/01/01"[PubDate])'
        TERM1 = f'("{Words[0]}"[Title/Abstract]) AND ({BegD}:{EndD}[Date - Publication])'
        # TERM1 = '('+'Words[0]'+'[Title/Abstract]) AND ' + '"' + BegD + '"' +'[PubDate] : ' + '"' +EndD + '"' +'[PubDate])'
        # первая версия: TERM1 = '('+'Words[0]'+'[ALL FIELDS]) AND ' + '"' + BegD + '"' +'[PubDate] : ' + '"' +EndD + '"' +'[PubDate])'
        # term = f'("{query}"[Title/Abstract]) AND ({start_date}:{end_date}[Date - Publication])'
#        print(TERM)
    elif len(Words) == 2:
        TERM1 = '(('+ Words[0] + '[ALL FIELDS]) AND ' + Words[1]+'[ALL FIELDS]) AND ' + '"' + BegD + '"' +'[PubDate] : ' + '"' +EndD + '"' +'[PubDate])'
#        print(TERM)
    elif len(Words) == 3: 
        TERM1 = '((('+ Words[0] + '[ALL FIELDS]) AND ' + Words[1]+'[ALL FIELDS]) AND ' + Words[2]+'[ALL FIELDS]) AND ' + '"' + BegD + '"' +'[PubDate] : ' + '"' +EndD + '"' +'[PubDate])'
#        print(TERM)
    return(TERM1)    
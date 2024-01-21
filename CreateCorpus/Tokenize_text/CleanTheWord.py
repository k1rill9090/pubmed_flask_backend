# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 02:24:09 2021

@author: PC
"""

def CleanWord(Word,LenWW):
    Delimiter = ''
    while LenWW > 0 and Word[0] in '/;.,!:"()[]': 
        Word = Word[1:LenWW]
        LenWW = LenWW - 1
    while LenWW > 0 and Word[-1] in '/;).,!:"()[]': 
        Word = Word[:-1]
        LenWW = LenWW - 1
        Delimiter = 'EndOfTerm'
    while LenWW > 0:
        if Word[-1] == ' ':
            Word = Word[:-1]
            LenWW = LenWW - 1
        else:
            break
    return(Word,Delimiter)

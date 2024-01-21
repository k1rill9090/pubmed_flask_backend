# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 02:22:32 2021

@author: PC
"""
from CreateCorpus import DB_PATH
import sys
import sqlite3
def read_Record_from_Corpus(idArt): # ,TypeStruct):
    try:
       conn = sqlite3.connect(DB_PATH)
       cursor = conn.cursor()
       sqlite_select_with_param = """SELECT Author,Title, Abstract,Year from ArticleStruct WHERE idArt = ?;"""
       cursor.execute(sqlite_select_with_param,(idArt,)) #, data_tuple)
       Res = cursor.fetchone()
       conn.commit()
       cursor.close()
       if Res != None: 
           return (Res[0],Res[1],Res[2],Res[3]) # return [Str, Year]
       else:
           return 'EndOfCorpus'
    except sqlite3.Error as error:
        print("Ошибка при работе с SQLite", error)
        raise
        # return('')
    finally:
        if conn:
            conn.close()
#            print("Соединение с SQLite закрыто")

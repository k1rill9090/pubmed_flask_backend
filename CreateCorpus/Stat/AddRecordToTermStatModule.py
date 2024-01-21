# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 13:02:16 2021

@author: PC
"""
from CreateCorpus import DB_PATH
# AddRecordToTermStatModule
import sqlite3
def AddRecordToTermStat(idTerm,StatNumber,Year,Term):
#    print('PRKey = ',PRKey,'  idTerm = ',idTerm,' StatNumber = ', StatNumber,'  Year = ',Year,'  Term = ', Term)
#    wait = input("PRESS ENTER TO CONTINUE.COUNT 88888888888")

    try:
       conn = sqlite3.connect(DB_PATH) # или :memory: чтобы сохранить в RAM
       cursor = conn.cursor()
       sqlite_insert_with_param = """INSERT INTO StatResult
                              (idTerm,StatNumber,Year,Term)
                              VALUES (?, ?, ?, ?);"""
       data_tuple = (idTerm,StatNumber,Year,Term)

#       print('SAVE 6666 - FirstTerm =',FirstTerm,' RestTerms =',RestTerms)
#       wait = input("SAVE - PRESS ENTER TO CONTINUE.")

       cursor.execute(sqlite_insert_with_param, data_tuple)
       conn.commit()
       cursor.close()

    except sqlite3.Error as error:
        print("Ошибка при работе с SQLite", error)
        raise
    finally:
        if conn:
            conn.close()
            print("Соединение с SQLite закрыто")
    
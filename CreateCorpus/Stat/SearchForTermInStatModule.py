# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 01:13:44 2021

@author: PC
"""
from CreateCorpus import DB_PATH
# import sys
import sqlite3
def SearchForTermIdInStat(Term1,Year1):
    try:
       conn = sqlite3.connect(DB_PATH) # или :memory: чтобы сохранить в RAM
       cursor = conn.cursor()
       sqlite_select_with_param = """SELECT Term, Year from StatResult WHERE Term = ? AND Year = ?;"""
       cursor.execute(sqlite_select_with_param,(Term1,Year1,))
       Res = cursor.fetchone() # Res = cursor.fetchall()
       conn.commit()
       cursor.close()
       
#       print('Res = ',Res)
#       wait = input("PRESS ENTER TO CONTINUE.COUNT 666666")
          
       if Res != None: 
           return (Res) # return [Str, Year]
       else:
#           return ('','NoSuchTerm')
           return ('NoSuchTerm')
    except sqlite3.Error as error:
        print("Ошибка при работе с SQLite", error)
        # return('')
        raise
    finally:
        if conn:
            conn.close()
            print("Соединение с SQLite закрыто")
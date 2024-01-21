# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 01:13:44 2021

@author: PC
"""
from CreateCorpus import DB_PATH
# import sys
import sqlite3
def TermsCountStat(Term1,Year1):
  try:
    
#    print('Term1 = ',Term1,' Year1 = ',Year1)
#    wait = input("PRESS ENTER TO CONTINUE.COUNT 666666")
      
    
    conn = sqlite3.connect(DB_PATH) # или :memory: чтобы сохранить в RAM
    cursor = conn.cursor()
    sqlite_select_with_param = """SELECT COUNT(*) as count from Terms WHERE Term = ? AND Year = ?;"""
    cursor.execute(sqlite_select_with_param,(Term1,Year1,))
    Res = cursor.fetchone() # Res = cursor.fetchall()
    conn.commit()

#    print('Res = ',Res,' Res[0] = ',Res[0])
#    wait = input("PRESS ENTER TO CONTINUE.COUNT 666666")

    cursor.close()
    return (Res[0]) # return [Str, Year]
  except sqlite3.Error as error:
    print("Ошибка при работе с SQLite", error)
    # return('')
    raise
  finally:
    if conn:
      conn.close()
      print("Соединение с SQLite закрыто")       


# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 19:30:04 2021

@author: PC
"""
from CreateCorpus import DB_PATH
# import sys
import sqlite3
def ReadTermFromTermsTable(idTerm1):

#    if idTerm1 > 168:
#       print('ReadTermFromTermsTable, idTerm1 = ',idTerm1)
#       wait = input("PRESS ENTER TO CONTINUE.121212121212")
    
    try:
       conn = sqlite3.connect(DB_PATH) # или :memory: чтобы сохранить в RAM
       cursor = conn.cursor()
       sqlite_select_with_param = """SELECT Term,Year from Terms WHERE idTerm = ?;"""
       cursor.execute(sqlite_select_with_param,(idTerm1,))
#       Res = cursor.fetchall()
       Res = cursor.fetchone()

#    if idTerm1 > 168:
#       print('idTerm1 = ',idTerm1,'  Res = ',Res)
#       wait = input("PRESS ENTER TO CONTINUE.23232323232323")
       
#       print('111')

#       wait = input("PRESS ENTER TO CONTINUE.")
       
       conn.commit()
       cursor.close()
       
#       print('222 Res = ', Res)
       
       if Res != None: 
           print('333 - Res = ',Res)
           return (Res[0],Res[1]) # return [Str, Year]
       else:
           return ('EndOfTerms','')
       
       print('444')
        
    except sqlite3.Error as error:
        print("Ошибка при работе с SQLite", error)
        # return('')
        raise
    finally:
        if conn:
            conn.close()
            print("Соединение с SQLite закрыто")
            
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 03:00:42 2021

@author: PC
"""
from CreateCorpus import DB_PATH
import sqlite3
def save_TermToTheDict(Term,TermAddr,TypeStruct,LenTerm,idArt,Year):
    try:
       conn = sqlite3.connect(DB_PATH)
       cursor = conn.cursor()
       sqlite_insert_with_param = """INSERT OR IGNORE INTO Terms
                              (TermAddr, Term, TypeStruct, LenTerm, idArt,Year)
                              VALUES (?, ?, ?, ?, ?, ?);"""
    #    idTerm = idTerm + 1 
       data_tuple = (TermAddr, Term, TypeStruct, LenTerm, idArt, Year)
       cursor.execute(sqlite_insert_with_param, data_tuple)
       conn.commit()
       cursor.close()

    except sqlite3.Error as error:
        print("Ошибка при работе с SQLite", error)
        raise
    finally:
        if conn:
            conn.close()
#            print("Соединение с SQLite закрыто")

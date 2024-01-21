#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from CreateCorpus import DB_PATH
import sqlite3
def insert_varible_into_table(Author, Title, Abstract, Year):
    try:
#       print('idArt = ',idArt,' Author = ',Author, ' Title = ',Title, ' Abstract = ', ' Year = ',Year)
#       conn = sqlite3.connect("C:\Test2\PubMedArticles-7.db") # или :memory: чтобы сохранить в RAM
       conn = sqlite3.connect(DB_PATH) # или :memory: чтобы сохранить в RAM
#       print('222222222')
       cursor = conn.cursor()
#       print('3333333333333')
       sqlite_insert_with_param = """INSERT OR IGNORE INTO ArticleStruct
                              (Author, Title, Abstract, Year)
                              VALUES (?, ?, ?, ?);"""

#       print('4444444')
       data_tuple = (Author, Title, Abstract, Year)
#       print('555555')

       cursor.execute(sqlite_insert_with_param, data_tuple)
#       print('666666')

       conn.commit()
       cursor.close()

    except sqlite3.Error as error:
       #  print("Ошибка при работе с SQLite", error)
       raise
    finally:
        if conn:
            conn.close()
            print("Соединение с SQLite закрыто")

# insert_varible_into_table(1, 'Author', 'Title', 'Abstract', 'Year')



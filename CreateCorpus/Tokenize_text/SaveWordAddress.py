import sqlite3
from CreateCorpus import DB_PATH
def save_WordAddr(WordAddr,LenW,idArt):
    try:
       conn = sqlite3.connect(DB_PATH)
       cursor = conn.cursor()
       sqlite_insert_with_param = """INSERT OR IGNORE INTO MarkWords
                              (WordAddr, LenW, idArt)
                              VALUES ( ?, ?, ?);"""
       data_tuple = (WordAddr,LenW, idArt)
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

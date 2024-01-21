import sqlite3
from CreateCorpus import DB_PATH

def clear_db():
    print("Запуск процедуры очистки БД")
    # обратить внимание на путь, если проблемы с соединением к бд
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # querry = "SELECT name FROM sqlite_master WHERE type='table';"

    remove_articles = "DELETE FROM ArticleStruct"
    cursor.execute(remove_articles)
    # очистить значение автоинкремента (по умолчанию не сбрасывается после очистки таблицы)
    cursor.execute('DELETE FROM sqlite_sequence WHERE name = "ArticleStruct"')
    conn.commit()
    cursor.execute("SELECT * FROM ArticleStruct")
    print("записи в бд ArticleStruct:", cursor.fetchall(), sep="\n")

    remove_markWords = "DELETE FROM MarkWords"
    cursor.execute(remove_markWords)
    cursor.execute('DELETE FROM sqlite_sequence WHERE name = "MarkWords"')
    conn.commit()
    cursor.execute("SELECT * FROM MarkWords")
    print("записи в бд MarkWords:", cursor.fetchall(), sep="\n")

    remove_statResult = "DELETE FROM StatResult"
    cursor.execute(remove_statResult)
    cursor.execute('DELETE FROM sqlite_sequence WHERE name = "StatResult"')
    conn.commit()
    cursor.execute("SELECT * FROM StatResult")
    print("записи в бд StatResult:", cursor.fetchall(), sep="\n")

    # remove_statResult2 = "DELETE FROM StatResult2"
    # cursor.execute(remove_statResult2)
    # conn.commit()
    # cursor.execute("SELECT * FROM StatResult2")
    # print("записи в бд StatResult2:", cursor.fetchall(), sep="\n")

    # remove_termStat = "DELETE FROM TermStat"
    # cursor.execute(remove_termStat)
    # conn.commit()
    # cursor.execute("SELECT * FROM TermStat")
    # print("записи в бд TermStat:", cursor.fetchall(), sep="\n")

    remove_terms = "DELETE FROM Terms"
    cursor.execute(remove_terms)
    cursor.execute('DELETE FROM sqlite_sequence WHERE name = "Terms"')
    conn.commit()
    cursor.execute("SELECT * FROM Terms")
    print("записи в бд Terms:", cursor.fetchall(), sep="\n")

    conn.close()

def check_db():
    print("проверка записей в БД в таблице ArticleStruct")
    # обратить внимание на путь, если проблемы с соединением к бд
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # querry = "SELECT name FROM sqlite_master WHERE type='table';"

    # посмотреть названия таблиц в БД
    # cursor.execute(querry)

    cursor.execute("SELECT * FROM ArticleStruct")
    answer = cursor.fetchall()

    conn.close()
    return answer

def check_terms():
    '''проверка таблицы Terms на наличие записей'''
    print("проверка записей в таблице Terms")
    # обратить внимание на путь, если проблемы с соединением к бд
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # querry = "SELECT name FROM sqlite_master WHERE type='table';"

    # посмотреть названия таблиц в БД
    # cursor.execute(querry)

    cursor.execute("SELECT * FROM Terms")
    answer = cursor.fetchall()

    conn.close()
    return answer

def check_article_exist(title='tittle'):
    '''проверка на наличие записи в таблице ArticleStruct, если ответ True, значит такой записи нет, можно записывать в БД'''
    # обратить внимание на путь, если проблемы с соединением к бд
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # querry = "SELECT name FROM sqlite_master WHERE type='table';"

    cursor.execute(f'''SELECT title FROM ArticleStruct WHERE title = "{title}"''')
    answer = cursor.fetchall()
    conn.close()
    # print("записи в бд:", answer, sep="\n")
    if len(answer) == 0:
        return True
    else:
        return False

def check_markWords_exist(wordAddr=0, lenW=3, idArt=1):
    '''Функция которая проверяет наличие записи в таблице MarkWords. Если такой записи нет, то выводится True, можно записывать в таблицу'''

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # querry = "SELECT name FROM sqlite_master WHERE type='table';"

    cursor.execute(f'''SELECT idArt, WordAddr, LenW FROM MarkWords WHERE idArt={idArt} AND WordAddr = {wordAddr} and LenW = {lenW}''')
    answer = cursor.fetchall()
    conn.close()
    # print("записи в бд:", answer, sep="\n")
    if len(answer) == 0:
        return True
    else:
        return False

def check_terms_exist(term, termAddr, typeStruct, lenTerm, idArt, year):
    '''Функция которая проверяет наличие записи в таблице Terms. Если такой записи нет, то выводится True, можно записывать в таблицу'''

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # querry = "SELECT name FROM sqlite_master WHERE type='table';"

    cursor.execute(f'''SELECT Term, TermAddr, TypeStruct, LenTerm, idArt, Year FROM Terms WHERE Term = "{term}" AND TermAddr = {termAddr} AND TypeStruct = {typeStruct} AND LenTerm = {lenTerm} AND idArt = {idArt} AND Year = {year}''')
    answer = cursor.fetchall()
    conn.close()
    # print("записи в бд:", answer, sep="\n")
    if len(answer) == 0:
        return True
    else:
        return False
    
def check_statResult_exist(idTerm, statNumber, year, term):
    '''Функция которая проверяет наличие записи в таблице StatResult. Если такой записи нет, то выводится True, можно записывать в таблицу'''

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # querry = "SELECT name FROM sqlite_master WHERE type='table';"

    cursor.execute(f'''SELECT idTerm, StatNumber, Year, Term FROM StatResult WHERE idTerm={idTerm} AND StatNumber = {statNumber} AND Year = {year} AND Term = "{term}"''')
    answer = cursor.fetchall()
    conn.close()
    # print("записи в бд:", answer, sep="\n")
    if len(answer) == 0:
        return True
    else:
        return False

def get_idArt():
    '''Функция, которая находит значение idArt первой записи в таблице ArticleStruct'''
    # обратить внимание на путь, если проблемы с соединением к бд
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # querry = "SELECT name FROM sqlite_master WHERE type='table';"

    cursor.execute(f"SELECT idArt FROM ArticleStruct ORDER BY idArt ASC LIMIT 1")
    answer = cursor.fetchall()
    conn.close()
    return answer[0][0]

def get_idTerm():
    '''Функция, которая находит значение idTerm последней записи в таблице Terms'''
    # обратить внимание на путь, если проблемы с соединением к бд
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # querry = "SELECT name FROM sqlite_master WHERE type='table';"

    cursor.execute(f"SELECT idArt FROM ArticleStruct ORDER BY idArt DESC LIMIT 1")
    answer = cursor.fetchall()
    conn.close()
    return answer[0][0]

if __name__ == "__main__":
    # print(check_db())
    print(check_terms())
    # a = check_article_exist()
    # print(get_idArt())
    # print(get_idTerm())
    # clear_db()
    # print(check_markWords_exist())
    # print(check_statResult_exist(2, 2, 2023, 'increase'))
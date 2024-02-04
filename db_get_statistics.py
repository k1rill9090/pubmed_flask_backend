import sqlite3
import json
from CreateCorpus import DB_PATH


def get_statistics(limit, offset, term, year):
    print("получение терминов из StatResult")

    # если LIMIT = None, то по умолчанию ставить 10 (формировать строку для запроса)
    if limit == None:
        limit = 100
        limit_str = f'LIMIT {limit}'
    else:
        limit_str = f'LIMIT {limit}'

    # если OFFSET = None, то формировать строку без  OFFSET
    if offset == None:
        offset = 0
        offset_str = f'OFFSET {offset}'
    else:
        offset_str = f'OFFSET {offset}'

    if term != None or year != None:
        where_str = "WHERE"
        if term != None and year != None:
            and_str = "AND"
        else:
            and_str = ""
    else:
        where_str = ""
        and_str = ""
    # если term = None, то формировать строку без  отбора по полю StatResult.Term
    if term == None:
        term_str = ''
    else:
        term_str = f"Term like '%{term}%' "

    # если year = None, то формировать строку без  отбора по полю StatResult.Term
    if year == None:
        year_str = ''
    else:
        year_str = f"year = {year}"

    # обратить внимание на путь, если проблемы с соединением к бд
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # запрос для вывода общего кол-ва записей
    cursor.execute(f"SELECT count(*) from StatResult")
    total_count = cursor.fetchall()[0][0]

    # основной запрос на получение данных
    cursor.execute(f"SELECT Term, StatNumber, Year FROM StatResult {where_str} {term_str} {and_str} {year_str} ORDER by Year ASC, StatNumber DESC {limit_str} {offset_str}")

    ans = cursor.fetchall()

    
    # print(len(ans))
    ans_dct = {
        "meta": {
            'limit': limit,
            'offset': offset,
            'total_count': total_count
        },
        "data": list()
    }

    for elem in ans:
        ans_dct['data'].append(
            {
            'termName': elem[0],
            'numOfAppearance': elem[1],
            'year': elem[2]
            }
        ) 

    
    ans_json = json.dumps(ans_dct, ensure_ascii=False)
    conn.close()
    return ans_json
if __name__ == "__main__":
    get_statistics()
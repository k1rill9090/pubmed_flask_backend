import sqlite3
import json
from CreateCorpus import DB_PATH


def get_terms(limit, offset):
    print("получение терминов из Terms")

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

    # обратить внимание на путь, если проблемы с соединением к бд
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    cursor.execute(f"SELECT count(*) from Terms")
    count = cursor.fetchall()[0][0]

    cursor.execute(f"SELECT idTerm, TermAddr, Term, idArt FROM Terms {limit_str} {offset_str}")

    ans = cursor.fetchall()

    
    # print(len(ans))
    ans_dct = {
        "meta": {
            'limit': limit,
            'offset': offset,
            'total_count': count
        },
        "data": list()
    }

    for elem in ans:
        ans_dct['data'].append(
            {
            'id': elem[0],
            'termAddr': elem[1],
            'termName': elem[2],
            'idArticle': elem[3]
            }
        ) 

    
    ans_json = json.dumps(ans_dct, ensure_ascii=False)
    conn.close()
    return ans_json
if __name__ == "__main__":
    get_terms()
import sqlite3
import json
from CreateCorpus import DB_PATH


def get_art(limit, offset):
    print("получение статей из ArticleStruct")

    # если LIMIT = None, то по умолчанию ставить 10 (формировать строку для запроса)
    if limit == None:
        limit = 10
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

    cursor.execute(f"SELECT count(*) from ArticleStruct")
    count = cursor.fetchall()[0][0]

    cursor.execute(f"SELECT idArt, title, abstract FROM ArticleStruct {limit_str} {offset_str}")

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
            'title': elem[1],
            'article': elem[2]
            }
        ) 

    
    ans_json = json.dumps(ans_dct, ensure_ascii=False)
    conn.close()
    return ans_json
if __name__ == "__main__":
    get_art()
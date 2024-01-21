import sqlite3
import json
from CreateCorpus import DB_PATH


def get_art(limit, offset):
    print("получение статей из ArticleStruct")
    # обратить внимание на путь, если проблемы с соединением к бд
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    cursor.execute(f"SELECT count(*) from ArticleStruct")
    count = cursor.fetchall()[0][0]

    cursor.execute(f"SELECT idArt, title, abstract FROM ArticleStruct {f'LIMIT {limit}' if limit != None else ''} {f'OFFSET {offset}' if offset != None else ''}")

    ans = cursor.fetchall()
    # print(len(ans))
    ans_dct = {
        "meta": {
            'limit': limit,
            'offset': offset,
            'count': count
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
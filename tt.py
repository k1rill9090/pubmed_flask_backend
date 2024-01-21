from Bio import Entrez
import xml.etree.ElementTree as ET
from Bio import Medline

def search_pubmed(query):
    Entrez.email = 'ol-zolot@yandex.ru'  # Укажите ваш адрес электронной почты

    handle = Entrez.esearch(db='pubmed', term=query, retmax=5)  # Укажите максимальное количество статей (в данном случае 10)
    record = Entrez.read(handle)
    handle.close()

    id_list = record['IdList']
    return id_list

def fetch_fulltext(id_list):
    print("список id: ", id_list)
    handle = Entrez.efetch(db='pubmed', id=id_list, rettype='medline', retmode='text')

    ret = Medline.parse(handle)
    return ret

# Пример использования

query = f'("covid-19 breast cancer"[Title/Abstract]) AND (2015/01/01:2020/02/01[Date - Publication])'
id_list = search_pubmed(query)
articles = fetch_fulltext(','.join(id_list))
for i in articles:
    try:
        print(i['AB'])
    except KeyError as exp:
        print('ошбочка', type(exp))
        if type(exp) == KeyError:
            print('2121')
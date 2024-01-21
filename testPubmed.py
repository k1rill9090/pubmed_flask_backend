from Bio import Entrez
import xml.etree.ElementTree as ET

def search_pubmed(query):
    Entrez.email = 'ol-zolot@yandex.ru'  # Укажите ваш адрес электронной почты

    handle = Entrez.esearch(db='pubmed', term=query, retmax=1)  # Укажите максимальное количество статей (в данном случае 10)
    record = Entrez.read(handle)
    handle.close()

    id_list = record['IdList']
    return id_list

def fetch_fulltext(id_list):
    print("список id: ", id_list)
    handle = Entrez.efetch(db='pubmed', id=id_list, rettype='XML', retmode='xml')
    results = handle.read()
    handle.close()

    return results

# Пример использования

query = f'("covid"[Title/Abstract]) AND (1800-01-01:2020/02/01[Date - Publication])'
id_list = search_pubmed(query)
articles = fetch_fulltext(','.join(id_list)).decode('utf-8')

print(articles)

articles_data = []
root_node = ET.ElementTree(ET.fromstring(articles)).getroot()
# print(root_node)
# формируем массив, содержащий словари с данными по каждой загруженной статье
for articles in root_node:
    data = dict()
    data['title'] = articles.find('MedlineCitation/Article/ArticleTitle').text
    # print("\n=======\n", data['title'])

    abstracts = []
    # поиск всех абстрактов (в некоторых статьях их может быть несколько)
    for tag in articles.findall('MedlineCitation/Article/Abstract/AbstractText'):
        # получить значение атрибута у тега
        label = tag.get('Label')
        # print(label)
        text = tag.text
        if label != None:
            abstracts.append(f'\n{label}\n{text}')
        else:
            abstracts.append(text)
    abstracts_text = '\n'.join(abstracts)
    data['abstract'] = abstracts_text
    # print(data['abstract'])

    date_year = articles.find('MedlineCitation/Article/ArticleDate/Year').text
    date_month = articles.find('MedlineCitation/Article/ArticleDate/Month').text
    date_day = articles.find('MedlineCitation/Article/ArticleDate/Day').text
    date = f'{date_day}.{date_month}.{date_year}'
    data['date'] = date
    # print(data['date'])


    last_name_arr = []
    # поиск всех авторов
    for tag in articles.findall('MedlineCitation/Article/AuthorList/Author/LastName'):
        # формирование массива из всех фамилий авторов
        last_name_arr.append(tag.text)
        # print(tag.text)
    initials_arr = []
    # поиск всех авторов
    for tag in articles.findall('MedlineCitation/Article/AuthorList/Author/Initials'):
        # формирование массива из всех инициалов авторов
        initials_arr.append(tag.text)
        # print(tag.text)
    authors = []
    # смешиваем фамилии и инициалы в один элемент массива
    for i in range(0, len(last_name_arr)-1):
        authors.append(f'{last_name_arr[i]} {initials_arr[i]}')

    authors_text = ', '.join(authors)
    data['authors'] = authors_text
    # print('+++', data['authors'])
    articles_data.append(data)

print(articles_data)

#!/usr/bin/env python
# coding: utf-8

# In[14]:


from Bio import Entrez
from Bio import Medline
from CreateCorpus import RecordToDbModule
from CreateCorpus import FormStrModule
# from . import check_delete_db
# import os
# import sys
# sys.path.append(".")
# from CreateCorpus import RecordToDbModule
import parse_pubmed_xml
import time



def create_corpus(count: int, search_line: str, start_date: str, end_date: str):
    '''обернул код в этом модуле в процедуру, чтобы его можно было выполнять из другого модуля(вызывается в app.py).
    Принимает на вход параметры, полученные в app.py из html формы'''

    print('\n==========================================\nЗапуск модуля Create_corpus')

    Entrez.email = 'ol-zolot@yandex.ru'

    # TERM='breast cancer' # Запрос к PubMed
    # TERM='(biomarker) AND (personalized medicine) AND (breast cancer)' # Запрос к PubMed
    # TERM = '(("coronavirus"[MeSH Terms] OR "coronavirus"[All Fields]) AND ("COVID-19"[All Fields] OR "COVID-19"[MeSH Terms] OR "COVID-19 Vaccines"[All Fields] OR "COVID-19 Vaccines"[MeSH Terms] OR "COVID-19 serotherapy"[All Fields] OR "COVID-19 Nucleic Acid Testing"[All Fields] OR "covid-19 nucleic acid testing"[MeSH Terms] OR "COVID-19 Serological Testing"[All Fields] OR "covid-19 serological testing"[MeSH Terms] OR "COVID-19 Testing"[All Fields] OR "covid-19 testing"[MeSH Terms] OR "SARS-CoV-2"[All Fields] OR "sars-cov-2"[MeSH Terms] OR "Severe Acute Respiratory Syndrome Coronavirus 2"[All Fields] OR "NCOV"[All Fields] OR "2019 NCOV"[All Fields] OR (("coronavirus"[MeSH Terms] OR "coronavirus"[All Fields] OR "COV"[All Fields]) AND 2019/11/01[PubDate] : 3000/12/31[PubDate]))) AND ("2020/10/01"[PubDate] : "2020/10/31"[PubDate])'
    # TERM = '(coronavirus) and (2020/10/01 : 2020/10/31)'
    # TERM = '((coronavirus) OR (coronavirus"[MeSH Terms] OR(coronavirus"[All Fields])) AND (2020/10/01 : 2020/10/31))'

    # Попытка соединиться в PMC
    # TERM = '("coronavirus"[MeSH Terms] OR "coronavirus"[All Fields]) AND ("2020/10/01"[PubDate] : "2020/10/31"[PubDate]) '
    # TERM = '("coronavirus"[MeSH Terms] OR "coronavirus"[All Fields]) and ("2021/12/17"[Date - Publication] : "2021/12/24"[Date - Publication])'
    # TERM = '("covid"[MeSH Terms] OR "covid"[All Fields]) and ("insomnia"[MeSH Terms] OR "insomnia"[All Fields]) and ("2020/01/01"[Date - Publication] : "2022/03/12"[Date - Publication])'

    # 42 записи только вместо 476
    # TERM = '(("breast neoplasms"[MeSH Terms] OR ("breast"[All Fields] AND "neoplasms"[All Fields]) OR "breast neoplasms"[All Fields] OR ("breast"[All Fields] AND "cancer"[All Fields]) OR "breast cancer"[All Fields]) AND ("calcinosis"[MeSH Terms] OR "calcinosis"[All Fields] OR "microcalcification"[All Fields]) AND "cluster"[All Fields]) AND ("2012/01/01"[Date - Publication] : "2021/12/31"[Date - Publication])'

    # 16 записей
    # TERM = '("breast" AND "cancer" AND "microcalcification" AND "cluster") AND ("2012/01/01"[Print Publication Date] : "2021/31/12"[Print Publication Date])'

    # 16 записей
    # TERM = '((((cancer) AND breast) AND calcification) AND cluster) AND ("2012/01/01"[Print Publication Date] : "2021/31/12"[Print Publication Date])'

    # PMC - 754 записси
    # TERM = '((((cancer) AND breast) AND calcification) AND cluster) AND ("2012/01/01"[PrintPubDate] : "2021/31/12"[PrintPubDate])'

    # PMC -  записси
    # TERM = '((((cancer) AND breast) AND calcification) AND cluster) AND ("2012/01/01"[Date - Publication] : "2021/31/12"[Date - Publication])'

    # Для тестирования
    # TERM = '(((((antigen) AND cancer) AND breast) AND calcification) AND cluster) AND ("2020/01/01"[PrintPubDate] : "2021/31/12"[PrintPubDate])'

    # PMC - 647 записси
    # TERM = '((((cancer) AND breast) AND calcification) AND cluster) AND ("2012/01/01"[PrintPubDate] : "2021/31/12"[PrintPubDate])'

    # PMC - 938 статей
    # TERM = '((((cancer) AND breast) AND calcification) AND cluster) AND ("2004/01/01"[PrintPubDate] : "2021/31/12"[PrintPubDate])'

    # PMC - 6044 статей
    # TERM = '(((cancer) AND breast) AND calcification) AND ("2003/01/01"[PrintPubDate] : "2021/31/12"[PrintPubDate])'

    # PMC - 10001 статей - Максимальное количество = 10 000 - может ошибка с датами...
    #TERM = 'calcification AND cancer AND breast AND ("2003/01/01"[PubDate] : "2021/31/12"[PubDate])'

    # PMC - 418 статей - Максимальное количество = 10 000 - может ошибка с датами...
    # TERM = '(((calcification[Abstract]) AND cancer) AND breast) AND ("2003/01/01"[PrintPubDate] : "2021/31/12"[PrintPubDate])'

    # PMC -  841 статей - Максимальное количество = 10 000 - может ошибка с датами...
    TERM = FormStrModule.FormStr(search_line, start_date, end_date)

    MAX_COUNT = count # Ограничение количества выгруженных статей
    AU_Arr=''; TI=''; AB=''; DP=''

    try: h=Entrez.esearch(db='pubmed', retmax=MAX_COUNT, term=TERM)
    except: time.sleep(5);h=Entrez.esearch(db='pubmed', retmax=MAX_COUNT, term=TERM)
    result = Entrez.read(h)
    ids = result['IdList']
    handle = Entrez.efetch(db='pubmed', id=ids, rettype='XML', retmode='xml')
    result = handle.read().decode('utf-8')

    articles_data = parse_pubmed_xml.parse_xml(result)
    index = 1
    for re in articles_data:
        try:
            RecordToDbModule.insert_varible_into_table(index, re['authors'], re['title'], re['abstract'], re['date'])
        except Exception as exp:
            raise
        index += 1
        if index % 100:
            print(index)

if __name__ == "__main__":
    create_corpus()
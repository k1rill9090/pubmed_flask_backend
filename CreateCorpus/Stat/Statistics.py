# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 18:39:45 2021

@author: PC

Модуль для подсчета статистики
"""

from CreateCorpus.Stat import ReadTermFromTermsModule
from CreateCorpus.Stat import SearchForTermInStatModule
# import CalcStatModule
from CreateCorpus.Stat import AddRecordToTermStatModule
from CreateCorpus.Stat import TermsCountStatModule
import check_delete_db

def statistics ():

    if len(check_delete_db.check_terms()) == 0:
        print("Таблица с извлеченными терминами пустая. Выполнение данного модуля невозможно.")
        raise Exception("Таблица с извлеченными терминами пустая. Выполнение данного модуля невозможно.")
    '''Расчет статистики'''
    # indTerm - Индекс по основной таблице Term, полный обход, 1 раз
    # indRest - Индекс по таблице TermStat 
    # FirstTerm - адрес термина в таблице Term, записывается в TernmStat
    # Совпадающие значения в колонке FirstTerm - это количество вхождений термина 
    # FirstTerm в таблицу Terms - Статистика встречаемости
    indTerm = 0; Term = ''
    '''PRKey = 1 #задавать в таблице автоинкрементом'''
    while Term != 'EndOfTerms': # Конец записей в таблице, цикл по всем записям в БД
        Term = ''
        Year = ''
        indTerm = indTerm + 1   
        Res = ReadTermFromTermsModule.ReadTermFromTermsTable(indTerm) # Возврат термин и год
        if Res[0] == 'EndOfTerms':
            print('Term = ',Term,'   Res[0] = ',Res[0])    
            # wait = input("PRESS ENTER TO CONTINUE.111111")
            break
        else:
            Term = Res[0]
            Year = Res[1]
    
    #    wait = input("PRESS ENTER TO CONTINUE.222222")
        
        Res = SearchForTermInStatModule.SearchForTermIdInStat(Term,Year)
        
    #    print('Res = ',Res) # ,'   Res[1] = ',Res[1])    
    #    wait = input("PRESS ENTER TO CONTINUE.3333333")
        
    #    print('Res = ',Res)
    #    print('Term = ',Term,'  Year = ',Year)
    #    wait = input("PRESS ENTER TO CONTINUE.55555555555555")
        
        
        if Res == 'NoSuchTerm':

    #       print('Term = ',Term)    
        
            Res1 = TermsCountStatModule.TermsCountStat(Term,Year)

    #       print('Res1 = ',Res1)    
    #       wait = input("PRESS ENTER TO CONTINUE.777777777777")

            StatNumber = Res1 # Количество повторов
            if check_delete_db.check_statResult_exist(indTerm,StatNumber,Year,Term):
                AddRecordToTermStatModule.AddRecordToTermStat(indTerm,StatNumber,Year,Term)
            else: print('запись уже есть в таблице StatResult')
            '''PRKey = PRKey + 1 #задавать в таблице автоинкрементом'''
    #    if indTerm > 165:
    #       print('Term = ',Term)    
    #       wait = input("PRESS ENTER TO CONTINUE.222222")

            
    # Выбрать по очереди записи из Terms
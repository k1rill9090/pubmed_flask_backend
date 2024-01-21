#!/usr/bin/env python
# coding: utf-8

import nltk
import sys
sys.path.append('./CreateCorpus/Tokenize_text/')
from CreateCorpus.Tokenize_text import ReadRecordFromCorpus
from CreateCorpus.Tokenize_text import SaveTermToTheDict
from CreateCorpus.Tokenize_text import SaveWordAddress
from CreateCorpus.Tokenize_text import CleanTheWord
from CreateCorpus.Tokenize_text import ExtStopList
import check_delete_db


# перед выполнением этой функции добавить проверку на наличие записей в ArticleStruct,
# так как при  их отсутствии этот модуль выдаст ошибку
def Tokenize_text():
    '''функция, которая проводит токенизацию выгруженных текстов'''
    if len(check_delete_db.check_db()) == 0:
        print("Статьи не загружены. Выполнение данного модуля невозможно.")
        raise Exception("Статьи не загружены. Выполнение данного модуля невозможно.")

    global stopWords
    global currentYear

    def isThisATerm(WordStr):
        import nltk
        tokens = nltk.word_tokenize(WordStr)
        tagged = nltk.pos_tag(tokens)
        LenL = len(tagged)
        ind = 0
        while ind < LenL:
            if (('NN' in tagged[ind]) or ('FW' in tagged[ind]) or ('VBG' in tagged[ind])and(WordStr.lower() not in stopWords)):
                return(1)
            ind = ind + 1
        return(0)
    TypeStruct = 2 # Для отладки пока разбирается только аннотация - Abstract
    indArt = check_delete_db.get_idArt()-1 # Это idArt в таблице корпуса текстов
    '''idWord = 1 #добавлять в БД автоинкрементом'''
    S = ''
    '''idTerm = 0 #добавлять в БД автоинкрементом'''
    TypeStruct = 2
    currentYear = '2023'
    stopWords = ExtStopList.ExtStopWords()
    nltk.download('punkt')
    nltk.download('averaged_perceptron_tagger')
    while S != 'EndOfCorpus': # Конец записей в таблице, цикл по всем записям в БД
        i2 = 0
        indArt = indArt + 1   
        print('S = ',S)
        S = ReadRecordFromCorpus.read_Record_from_Corpus(indArt) # ,TypeStruct) 
        if S == 'EndOfCorpus':
            break
        Author = ''; Title = ''; Abstract = ''; Year = '';
        try: Author = S[0] 
        except: Author = '';
        try: Title = S[1] 
        except: Title = '';
        try: Abstract = S[2] 
        except: Abstract = '';
        try: 
            Year = S[3]
            if len(Year) >= 4: 
                Year = Year[0:4]
                currentYear = Year
            else:
                Year = currentYear
        except: Year = currentYear;
        if TypeStruct == 1: S = Title;
        if TypeStruct == 2: S = Abstract;
        LenS = len(S)
        Word = ''
        Term = ''
        Addr = 0
        isTerm = 0 # если в термине есть существительное или герундий, то это термин
        while i2 < LenS: 
            if S[i2].isalpha() and i2 < LenS:
                '''idWord = idWord + 1 #Добавлять в БД автоинкрементом'''
                Word = S[i2:LenS].split()[0].strip()
                LenW = len(Word)
                if check_delete_db.check_markWords_exist(i2, LenW, indArt):
                    SaveWordAddress.save_WordAddr(i2,LenW,indArt)
                else: print('данная запись уже есть в таблице MarkWords')
                Word = Word.lower()
                Word = CleanTheWord.CleanWord(Word,LenW)
                try: Delimiter = Word[1] 
                except: Delimiter = '';
                try: Word = Word[0] 
                except: Word = '';
                if Term == '':
                    Addr = i2
                    LenTerm = 0
                stopW = 0
                IsATerm = isThisATerm(Word)
                if Word.lower() not in stopWords and IsATerm == 1:
                    Term = Term + Word + ' '
                    LenTerm = LenTerm + 1
                else:
                    stopW = 1
                i2 =i2+len(Word) # Позиция следующего символа за словом
                LenOfStrinTerm = len(Term)
                if i2 < (LenS) and S[i2] != ' ' and Term != '' and not S[i2].isalpha() or stopW == 1 or Delimiter == 'EndOfTerm':
                    if isTerm == 0: isTerm = isThisATerm(Term)
                    if isTerm == 1: 
                        if LenOfStrinTerm > 3:
                            # условие такое же как и для MarkWords, так как эти данные связаны
                            # Без MarkWords не будет записей в Terms
                            if check_delete_db.check_terms_exist(Term.strip().lower(),Addr,TypeStruct,LenTerm,indArt,Year): 
                                SaveTermToTheDict.save_TermToTheDict(Term.strip().lower(),Addr,TypeStruct,LenTerm,indArt,Year)
                            else: print('данная запись уже есть в таблице Terms')
                            '''idTerm = idTerm + 1 #добавлять в БД автоинкрементом'''
                            Term = ''
                            isTerm = 0
                            stopW = 0
                            
                if (i2 >= (LenS)) and Term != ''  or Delimiter == 'EndOfTerm' or (i2 < (LenS))  and Term != '' and S[i2] != ' ' and not S[i2].isalpha():   
                    if isTerm == 0: 
                        isTerm = isThisATerm(Term)
                    if isTerm == 1: 
                        if LenOfStrinTerm > 3:
                            if check_delete_db.check_terms_exist(Term.strip().lower(),Addr,TypeStruct,LenTerm,indArt,Year):
                                SaveTermToTheDict.save_TermToTheDict(Term.strip().lower(),Addr,TypeStruct,LenTerm,indArt,Year)
                            else: print('данная запись уже есть в таблице Terms')
                            '''idTerm = idTerm + 1 #добавлять в БД автоинкрементом'''
                            Term = ''
                            isTerm = 0
                            stopW = 0
                            
            i2 = i2 + 1     
    file = open("./CreateCorpus/Tokenize_text/Params.txt", "w")
    file.write(str(check_delete_db.get_idTerm()))
    file.close()





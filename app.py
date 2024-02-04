from flask import Flask, render_template, request, redirect
from flask_cors import CORS, cross_origin
from CreateCorpus.CreateCorpus import create_corpus
from CreateCorpus.Tokenize_text import TokenizeText
from CreateCorpus.Stat import Statistics
from datetime import datetime
import db_get_articles
import db_get_terms
import db_get_statistics
import check_delete_db
import re
import json

from datetime import datetime

def checkInt(param: str) -> bool:
    '''Функция, которая проверяет полученное значение на целое число'''
    try:
        num = float(param)
        if int(str(num).split('.')[1]) == 0:
            return True
        else:
            return False
    except:
        return False

def checkDate(param: str) -> bool:
    '''Функция, которая проверяет полученное значение на корректность даты'''
    try:
        datetime.strptime(param, '%Y-%m-%d')
        return True
    except:
        return False

def checkQuerryParams(args, limit, offset, route=""):
    '''Функция, которая проверят корректность ввода querry параметров для методов GET /articles, /terms и /statistics
    параметры:
    args - querry-параметры, которые указываются при обращении к роуту
    limit и offset вынесены в отдельные параметры для доп. проверки
    route - название роута, к которому идет обращение (для роута /statistics прописана доп. логика)'''
     # готовим список ошибок
    errors = []
    for i in args:
        # если вызывается метод GET /statistics, то для него пропускать параметр term
        if route == 'statistics':
            if i == 'term' or i == 'year':
                continue
        if (i == 'limit' or i == 'offset'):
            if not(checkInt(args[i]) and int(args[i]) >= 0):
                errors.append(f"Параметр {i} должен быть целым числом")

        else:
            errors.append(f"Неизвестный параметр {i}")
    # проверка на запись одного поля offset
    if offset != None and limit == None:
        errors.append(f"Указан параметр offset без параметра limit. Укажите параметр limit!")
    return errors

app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'
app.config['JSON_AS_ASCII'] = False

# загрузить статьи в БД
@app.route("/articles", methods = ['POST'])
@cross_origin()
def postArticles():
    # готовим список обязательных полей
    reqFields = ['num', 'words', 'min_date', 'max_date']
    # готовим список ошибок
    errList = []
    for i in request.json:
        if i != 'data':
            errList.append(f'неизвестное поле {i}')
        else:
            for j in request.json[i]:
                if j not in reqFields:
                    errList.append(f'неизвестное поле в объекте {i}: {j}')

    for i in reqFields:
        if request.json['data'].get(i, None) == None:
            errList.append(f'В объекте data поле {i} обязательно для заполнения!')
        elif i == 'num' and not(checkInt(request.json['data'][i]) and int(request.json['data'][i]) > 0 and int(request.json['data'][i]) <= 100):
            errList.append(f'В объекте data поле {i} должно быть целым числом от 1 до 100')
        elif i == 'words' and not(request.json['data'][i] != '' and re.search(r'[^a-zA-Z\s\d-]+?', request.json['data'][i]) == None ):
            errList.append(f'В объекте data для поля {i} допускаются слова, разделенные пробелом и содержащие английские символы, цифры и тире')
        elif i == 'min_date' and not(request.json['data'][i] != '' and checkDate(request.json['data'][i]) and datetime.strptime(request.json['data'][i], '%Y-%m-%d') >= datetime.strptime('1800-01-01', '%Y-%m-%d')):
            errList.append(f'В объекте data поле {i} должно быть формата ГГГГ-ММ-ДД и не меньше 01.01.1800')
        elif i == 'max_date' and not(request.json['data'][i] != '' and checkDate(request.json['data'][i]) and datetime.strptime(request.json['data'][i], '%Y-%m-%d') >= datetime.strptime('1800-01-01', '%Y-%m-%d') ):
            errList.append(f'В объекте data поле {i} должно быть формата ГГГГ-ММ-ДД не меньше 01.01.1800')
        
    # доп проверка на дату
    if checkDate(request.json['data']['min_date']):
        if not(datetime.strptime(request.json['data']['min_date'], '%Y-%m-%d') <= datetime.strptime(request.json['data']['max_date'], '%Y-%m-%d') ):
            errList.append(f'В объекте data поле max_date должно быть не мешьше поля min_date')

    if len(errList) > 0:
        return json.dumps({ 'errors': errList }, ensure_ascii=False), 400, {'content-type':'application/json'}
    
    print(request.json)
    # распарсить json по переменным
    count_articles = request.json['data']['num']
    search = request.json['data']['words']
    start_date = datetime.strptime(request.json['data']['min_date'], '%Y-%m-%d').strftime('%Y/%m/%d')
    end_date = datetime.strptime(request.json['data']['max_date'], '%Y-%m-%d').strftime('%Y/%m/%d')
    # print(count_articles, search, start_date, end_date, sep='\n')

    print('выполнение модуля create_corpus')
    try:
        create_corpus(count_articles, search, start_date, end_date)
    except Exception as exp:
        errList.append(str(exp.args[0]))
    if len(errList) > 0:
        return json.dumps({ 'typeError': '500. Internal server error','errors': errList }, ensure_ascii=False), 500, {'content-type':'application/json'}
    
        # if type(exp) == KeyError:
        #     return json.dumps({'error': 'Записи до данному запросу не найдены'}, ensure_ascii=False), 500, {'content-type':'application/json'}
        # print("=-=-=-", exp.args)
        return json.dumps({'error':'Не удалось загрузить тексты из Pubmed', 'textErr': str(exp)}, ensure_ascii=False), 500, {'content-type':'application/json'}
    print('модуль create_corpus выполнен')
    return json.dumps({'status': 'ok', 'message': 'Тексты загружены'}, ensure_ascii=False,), 201, {'content-type':'application/json'}

# получить список загруженных статей
@app.route("/articles", methods = ['GET'])
@cross_origin()
def getArticles():
    # отдельно записать нужные querry параметры
    args = request.args
    limit = args.get('limit', None)
    offset = args.get('offset', None)
   
    errList = checkQuerryParams(args, limit, offset)
    # Если список ошибок непустой, то выдаем 400 ошибку и отправляем в response этот список 
    if len(errList) > 0:
        return json.dumps({'errors': errList}, ensure_ascii=False), 400, {'content-type':'application/json'}
    # если все ок, то возвращаем ответ с данными из бд
    print('формирование json Ответа')
    ans = db_get_articles.get_art(limit=limit, offset=offset)
    print('json готов')
    return ans, 200, {'content-type':'application/json'}

# Извлечение терминов из статей
@app.route('/terms', methods = ['POST'])
def postTerms():
    errList = []
    print('\n=====================================\nЗапуск модуля Tokenize_text')
    try:
        TokenizeText.Tokenize_text()
    except Exception as exp:
        errList.append(str(exp.args))
    if len(errList) > 0:
        return json.dumps({ 'typeError': '500. Internal server error','errors': errList }, ensure_ascii=False), 500, {'content-type':'application/json'}
        # raise
        # print(exp)
        # return json.dumps({'error': str(exp.args[0])}, ensure_ascii=False), 500, {'content-type':'application/json'}
    return json.dumps({'status': 'ok', 'message': 'Выделение терминов выполнено успешно'}, ensure_ascii=False,), 201, {'content-type':'application/json'}

# Получение списка терминов
@app.route('/terms', methods = ['GET'])
def getTerms():
    # отдельно записать нужные querry параметры
    args = request.args
    limit = args.get('limit', None)
    offset = args.get('offset', None)
   
    errList = checkQuerryParams(args, limit, offset)
    # Если список ошибок непустой, то выдаем 400 ошибку и отправляем в response этот список 
    if len(errList) > 0:
        return json.dumps({'errors': errList}, ensure_ascii=False), 400, {'content-type':'application/json'}
    # если все ок, то возвращаем ответ с данными из бд
    print('формирование json Ответа')
    ans = db_get_terms.get_terms(limit=limit, offset=offset)
    print('json готов')
    return ans, 200, {'content-type':'application/json'}

# Расчет статистики
@app.route('/statistics', methods = ['POST'])
def postStatistics():
    errList = []
    print('\n=====================================\nЗапуск модуля Statistics')
    try:
        Statistics.statistics()
    except Exception as exp:
        errList.append(str(exp.args))
    if len(errList) > 0:
        return json.dumps({ 'typeError': '500. Internal server error','errors': errList }, ensure_ascii=False), 500, {'content-type':'application/json'}
        # print(exp)
        # return json.dumps({'error': str(exp)}, ensure_ascii=False), 500, {'content-type':'application/json'}
    return json.dumps({'status': 'ok', 'message': 'Расчет статистики по терминам выполнен успешно!'}, ensure_ascii=False,), 201, {'content-type':'application/json'}

# Получить статистику по терминам
@app.route('/statistics', methods = ['GET'])
def getStatistics():
    # отдельно записать нужные querry параметры
    args = request.args
    limit = args.get('limit', None)
    offset = args.get('offset', None)
    term = args.get('term', None)
    year = args.get('year', None)
   
    errList = checkQuerryParams(args, limit, offset, route='statistics')
    # Если список ошибок непустой, то выдаем 400 ошибку и отправляем в response этот список 
    if len(errList) > 0:
        return json.dumps({'errors': errList}, ensure_ascii=False), 400, {'content-type':'application/json'}
    # если все ок, то возвращаем ответ с данными из бд
    print('формирование json Ответа')
    ans = db_get_statistics.get_statistics(limit=limit, offset=offset, term=term, year=year)
    print('json готов')
    return ans, 200, {'content-type':'application/json'}

# Очистить базу данных
@app.route('/clearPubmedArticles', methods = ['DELETE'])
def delPubmedArticles():
    try:
        check_delete_db.clear_db()
    except Exception as exp:
        return json.dumps({'error': str(exp)}, ensure_ascii=False), 500, {'content-type':'application/json'}

    print('json готов')
    return json.dumps({'message': 'Все таблицы очищены'}, ensure_ascii=False), 200, {'content-type':'application/json'}
   


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')

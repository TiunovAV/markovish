# dependencies are here
from datetime import datetime
from random import randint
import pandas as pd


# functions are here
def txt_seq_read(file, start='1', stop='2', include=True, printres=True):
    """func gets txt file with a sequence of acts
    convert to 2D array by stop and start symbols
    prints summary
    uppercase ignored
    args:
    start - leading symbol of a sequence (default '1')
    stop - ending symbol of a sequence (default '2')
    include - include or not start and stop to analysis (default True)
    printres - print summary or not (default True)"""
    with open(file) as doc:
        data = doc.read().split()
    data = [i.lower().strip() for i in data]  # convert to lowercase, remove spaces
    # check data for conditions
    if start not in data or data[0] != start:
        print('Check start symbol!')
    if stop not in data or data[-1] != stop:
        print('Check start symbol!')
    # create 2d arrray
    result = []
    box=[] # it worked well without this line but code inspector asked me to add it
    for i in range(len(data)):
        if data[i] == start:  # and i <= 1: # first case - we get 1st sequence
            box = []
            if include:
                box.append(data[i])
        elif data[i] == stop:
            if include:
                box.append(data[i])
            result.append(box)
        else:
            box.append(data[i])
    if printres:
        print(f'All data length = {len(data)}, includes {len(result)} sequences.')
    return result


def blanc_table(data, start='1', stop='2', order=True, include=False):
    """ creates blank table from 2D array, by default not include start and stop
    tries to order by occurrence in a sequence"""
    flatten = []
    for act in data:  # create flat array
        flatten.extend(act)
    labels = list(set(flatten))
    if not include:  # remove start and stop if needed
        res = []
        for i in labels:
            if i != start and i != stop:
                res.append(i)
        labels = res
    if order:  # order if chosen
        chart = dict.fromkeys(labels, 0)
        for i in labels:
            chart[i] = []
            for act in data:
                for j in range(len(act)):
                    if act[j] == i:
                        chart.get(act[j]).append(j)
            chart[i] = sum(chart[i]) / len(chart[i])
        chart = dict(sorted(chart.items(), key=lambda item: item[1]))
        labels = chart.keys()
    else:
        labels.sort()
    blanc = pd.DataFrame(0.0, index=labels, columns=labels)
    return blanc


def freq_table(blanc, data, flatten=False, saveres=True):
    """Creates markov transitions table to blanc tble from 
    2D array by default, or 1D array if flatten = True"""
    table = blanc.copy()
    if flatten:  # seems that this block is useless, but let it be here
        flatten = []
        for act in data:  # create flat array
            flatten.extend(act)
        data = flatten
        for start in table.index:  # getting events to compare
            for end in table.index:
                for i in range(len(data)):  # walk through data
                    if i + 1 < len(data):  # to prevent from out range error
                        # if we get data fit to headers - add 1 to a cell
                        if (data[i], data[i + 1]) == (start, end):
                            table.loc[start, end] = table.loc[start, end] + 1
    else:
        for start in table.index:
            for end in table.index:
                for act in data:
                    for i in range(len(act)):  # переберем данные
                        if i + 1 < len(act):  # если мы нашли в данных пару, которая соответствует заголовку столбца и строки
                            if (act[i], act[i + 1]) == (start, end):  # запишем единичку в значение данной ячейки
                                table.loc[start, end] = table.loc[start, end] + 1
    if saveres:
        fname = datetime.now().strftime('%d_%b_%Y_%H-%M') + '_transition_freq__table.xlsx'
        table.to_excel(fname)
    return table


def prob_table(blanc, trans, rd=2, saveres=True):
    """Receives blanc and raw transitions table (pd DataFrames)
    returns probabilities pd DataFrame
    args: round - nr of symbols after comma, default 2"""
    freq = blanc.copy()
    for start in trans.index:
        for end in trans.index:
            freq = freq.loc[:].astype(float)
            if trans.loc[start, :].sum() != 0:
                freq.loc[start, end] = round(trans.loc[start, end] / trans.loc[start, :].sum(), rd)
    if saveres:
        fname = datetime.now().strftime('%d_%b_%Y_%H-%M') + '_transitions_probability_table.xlsx'
        freq.to_excel(fname)
    return freq


def permutate(data, iters=9999):
    """gets 2D aray, returns array with each member copied by n = iters_of_permutations"""
    prognoz = len(data) * iters
    cnt = 0
    res = []
    print('Progress: ')
    for seq in data:
        for i in range(0, iters):
            replica = seq.copy()
            coin = randint(1, len(replica) - 3)
            replica[coin], replica[coin + 1] = replica[coin + 1], replica[coin]
            res.append(replica)
            print(round((cnt * 100) / prognoz, 2), end='\r')
            cnt += 1
    print('\r')
    print('Done!')
    return res


def get_seq_from_lens(file, sep=' - ', flatten=False):
    """func gets txt file with a sequence of acts seqs sep by 2\n after - time of event
    single string looks like:
    event - 05.66
    converts to 2D array 
    prints summary
    args:
    sep - specify separator
    flatten - turn 2D array to flat array
    """
    with open(file) as f:
        data = f.read().split('\n\n')  # режем по двойному разрыву строки
        data = [observ.split('\n') for observ in data]  # каждую строку превращаем в список
    result = []
    error_count = 0
    for observ in data:  # проходим по все наблюдениям из списка
        for line in observ:  # check data syntax
            condition = sep in line
            if not condition:
                print(f'Fix data - string {line} is invalid')
                error_count += 1
    if error_count > 0:
        raise BaseException('INVALID DATA')
    else:
        print("Data OK")

    for observ in data:
        boxes = ['begin']  # костыль
        # boxes = ['begin']
        for index, line in enumerate(observ):  # проходим по каждой строке наблюдения
            box = line.split(sep)  # получаем список из строки
            box = box[0]
            # строка готова, кладем её в корзиночку
            boxes.append(box)
            boxes.append('finish')  # еще один костыль
        result.append(boxes)  # когда корзиночка заполнена готовыми строками, кладем

    # уменьшаем уровень вложенности списков
    if flatten:
        flat = []
        for i in result:
            flat.extend(i)
        result = flatten

    return result

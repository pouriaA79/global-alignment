from itertools import combinations


def global_align(x, y, s_match, s_mismatch, s_gap):
    A = []
    for i in range(len(y) + 1):
        A.append([0] * (len(x) + 1))
    for i in range(len(y) + 1):
        A[i][0] = s_gap * i
    for i in range(len(x) + 1):
        A[0][i] = s_gap * i
    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            A[i][j] = max(
                A[i][j - 1] + s_gap,
                A[i - 1][j] + s_gap,
                A[i - 1][j - 1] + (s_match if (y[i - 1] == x[j - 1] and y[i - 1] != '-') else 0) + (
                    s_mismatch if (y[i - 1] != x[j - 1] and y[i - 1] != '-' and x[j - 1] != '-') else 0) + (
                    s_gap if (y[i - 1] == '-' or x[j - 1] == '-') else 0)
            )
    align_X = ""
    align_Y = ""
    i = len(x)
    j = len(y)
    while i > 0 or j > 0:
        current_score = A[j][i]
        if i > 0 and j > 0 and (
                ((x[i - 1] == y[j - 1] and y[j - 1] != '-') and current_score == A[j - 1][i - 1] + s_match) or
                ((y[j - 1] != x[i - 1] and y[j - 1] != '-' and x[i - 1] != '-') and current_score == A[j - 1][
                    i - 1] + s_mismatch) or
                ((y[j - 1] == '-' or x[i - 1] == '-') and current_score == A[j - 1][i - 1] + s_gap)
        ):
            align_X = x[i - 1] + align_X
            align_Y = y[j - 1] + align_Y
            i = i - 1
            j = j - 1
        elif i > 0 and (current_score == A[j][i - 1] + s_gap):
            align_X = x[i - 1] + align_X
            align_Y = "-" + align_Y
            i = i - 1
        else:
            align_X = "-" + align_X
            align_Y = y[j - 1] + align_Y
            j = j - 1
    return (align_X, align_Y, A[len(y)][len(x)])


def star_alignment(number_seq, seqs):
    tasks = tuple(combinations(zip(range(len(seqs)), seqs), 2))
    # print(seqs)
    matrix_s = [[0 for _ in range(number_s)] for _ in range(number_s)]
    matrix_a = [['' for _ in range(number_s)] for _ in range(number_s)]
    for i in range(len(tasks)):
        x, y, sco = global_align(tasks[i][0][1], tasks[i][1][1], 3, -1, -2)
        matrix_s[tasks[i][0][0]][tasks[i][1][0]] = sco
        matrix_a[tasks[i][0][0]][tasks[i][1][0]] = x + ' ' + y
        matrix_s[tasks[i][1][0]][tasks[i][0][0]] = sco
        matrix_a[tasks[i][1][0]][tasks[i][0][0]] = y + ' ' + x
    sum_of_matrix_score = []
    for i in range(number_s):
        sum_of_matrix_score.append(sum(matrix_s[i]))
    # print(matrix_s, matrix_a)
    center_index_r = sum_of_matrix_score.index(max(sum_of_matrix_score))
    new_list_center_row = list(set(matrix_s[center_index_r]))
    new_list_center_row.sort()
    order_of_a = []
    seen_zero = 0
    for j in reversed(range(len(new_list_center_row))):
        if matrix_s[center_index_r][
            matrix_s[center_index_r].index(new_list_center_row[j])] == 0 and seen_zero == 0:
            index_pos_list = [i for i in range(len(matrix_s[center_index_r])) if
                              matrix_s[center_index_r][i] == new_list_center_row[j]]
            for i in range(len(index_pos_list)):
                if index_pos_list[i] != center_index_r:
                    order_of_a.append(index_pos_list[i])
            seen_zero = 1
        else:
            index_pos_list = [i for i in range(len(matrix_s[center_index_r])) if
                              matrix_s[center_index_r][i] == new_list_center_row[j]]
            for i in range(len(index_pos_list)):
                order_of_a.append(index_pos_list[i])
    msa = []
    base = matrix_a[center_index_r][order_of_a[0]]
    base_align = base.split(' ')
    msa.append(base_align[0])
    msa.append(base_align[1])
    # print(order_of_a, 50, center_index_r)
    for t in range(number_seq - 2):
        # print(msa[0])
        next_item = seqs[order_of_a[t + 1]]
        # print(next_item)
        next_align, new_center, _ = global_align(next_item, msa[0], 3, -1, -2)
        msa.append(next_align)
        if new_center == msa[0]:
            continue
        else:
            # print(msa , 45)
            len_old_center = len(msa[0])
            # print(len_old_center)
            for j in range(t + 1):

                count = 0
                other = []
                for x in range(len_old_center):
                    if x < len(msa[j + 1]):
                        while True:
                            if x == 0 and new_center[count] == msa[0][x]:
                                other.append(msa[j + 1][x])
                                break
                            elif new_center[count + 1] == msa[0][x] and x != 0:
                                other.append(msa[j + 1][x])
                                count += 1
                                break
                            else:
                                other.append('-')
                                count += 1

                if count != len(new_center) - 1:
                    for yy in range(count + 1, len(new_center)):
                        # print(new_center[yy])
                        other.append(new_center[yy])
                msa[j + 1] = ''.join([str(elem) for elem in other])
            msa[0] = new_center

    # print(msa)
    correct_order_mas = []
    for i in range(number_seq):
        if i in order_of_a:
            correct_order_mas.append(msa[order_of_a.index(i) + 1])
        else:
            correct_order_mas.append(msa[0])
    return correct_order_mas


def calculate_score(msa):
    final_score = 0
    # print()
    for i in range(len(msa)):
        for j in range(i + 1, len(msa)):
            for t in range(len(msa[0])):
                # print(t)
                if msa[i][t] != msa[j][t] and msa[i][t] != '-' and \
                        msa[j][t] != '-':
                    final_score += -1
                elif msa[i][t] != msa[j][t] and msa[i][t] == '-':
                    final_score += -2
                elif msa[i][t] != msa[j][t] and msa[j][t] == '-':
                    final_score += -2
                elif msa[i][t] == msa[j][t] and msa[j][t] != '-':
                    final_score += 3
    return final_score


def block(old_msa, old_score):
    blocks = [["" for t in range(len(old_msa[0]))] for x in range(len(old_msa))]

    for i in range(len(old_msa[0])):
        words = []
        for j in range(len(old_msa)):
            words.append(old_msa[j][i])
        if len(set(words)) != 1:
            for j in range(len(old_msa)):
                blocks[j][i] = words[j]
    lines = []
    scores = []
    msas = []
    for i in range(len(old_msa[0]) - 1):
        if blocks[0][i] == "":
            lines.append(i)
    if len(lines) != 0:
        for i in range(len(lines) + 1):
            if i == 0 and lines[0] - 0 > 1:
                sqs = [["" for t in range(lines[0])] for x in range(len(old_msa))]
                count = 0
                for j in range(lines[0]):
                    for t in range(len(old_msa)):
                        if blocks[t][j] == '-':
                            sqs[t][count] = ""
                        else:
                            sqs[t][count] = blocks[t][j]
                    count += 1
                new_msas = []
                for t in range(len(old_msa)):
                    strr = ""
                    for j in range(len(sqs[t])):
                        if sqs[t][j] != "":
                            strr += sqs[t][j]
                    new_msas.append(strr)

                correct_msa = []
                new_align_ms = star_alignment(len(new_msas), new_msas)
                for t in range(len(old_msa)):
                    new_strr = new_align_ms[t]
                    old_strr = old_msa[t]
                    old_fraction = old_msa[t][0: lines[0]]
                    XXX = old_strr.replace(old_fraction, new_strr)
                    correct_msa.append(XXX)
                new_score = calculate_score(correct_msa)
                scores.append(new_score)
                msas.append(correct_msa)

            elif i > 0 and i < len(lines) and lines[i] - lines[i - 1] > 2:
                sqs = [["" for t in range(lines[i] - lines[i - 1] - 1)] for x in range(len(old_msa))]
                count = 0
                for j in range(lines[i - 1] + 1, lines[i]):
                    for t in range(len(old_msa)):
                        if blocks[t][j] == '-':
                            sqs[t][count] = ""
                        else:
                            sqs[t][count] = blocks[t][j]
                    count += 1

                new_msas = []
                for t in range(len(old_msa)):
                    strr = ""
                    for j in range(len(sqs[t])):
                        if sqs[t][j] != "":
                            strr += sqs[t][j]
                    new_msas.append(strr)

                correct_msa = []
                new_align_ms = star_alignment(len(new_msas), new_msas)

                for t in range(len(old_msa)):
                    new_strr = new_align_ms[t]
                    old_strr = old_msa[t]
                    old_fraction = old_msa[t][lines[i - 1] + 1: lines[i]]
                    XXX = old_strr.replace(old_fraction, new_strr)
                    correct_msa.append(XXX)

                new_score = calculate_score(correct_msa)
                scores.append(new_score)
                msas.append(correct_msa)

            elif i == len(lines) and len(old_msa[0]) - lines[i - 1] > 2:
                sqs = [["" for t in range(len(old_msa[0]) - lines[i - 1] - 1)] for x in range(len(old_msa))]
                count = 0

                for j in range(lines[i - 1] + 1, len(old_msa[0])):
                    for t in range(len(old_msa)):
                        if blocks[t][j] == '-':
                            sqs[t][count] = ""
                        else:
                            sqs[t][count] = blocks[t][j]
                    count += 1

                new_msas = []
                for t in range(len(old_msa)):
                    strr = ""
                    for j in range(len(sqs[t])):
                        if sqs[t][j] != "":
                            strr += sqs[t][j]
                    new_msas.append(strr)

                correct_msa = []
                new_align_ms = star_alignment(len(new_msas), new_msas)

                for t in range(len(old_msa)):
                    new_strr = ""
                    for j in range(len(old_msa[t])+100):
                        if j <= lines[i - 1]:
                            new_strr += old_msa[t][j]
                        else:
                            if j - lines[i - 1] - 1 < len(new_align_ms[t]):
                                new_strr += new_align_ms[t][j - lines[i - 1] - 1]

                    correct_msa.append(new_strr)

                new_score = calculate_score(correct_msa)
                scores.append(new_score)
                msas.append(correct_msa)

        if max(scores) > old_score:
            index = scores.index(max(scores))
            return scores[index], msas[index]
    return old_score, old_msa


if __name__ == "__main__":
    number_s = int(input())
    sqs = []
    for i in range(number_s):
        sq = input()
        sqs.append(sq)
    first_res = star_alignment(number_s, sqs)
    first_score = calculate_score(first_res)
    # print(first_score)
    new_score, new_msa = block(first_res, first_score)
    for j in range(100):
        if new_score > first_score:
            first_score = new_score
            first_res = new_msa
            new_score, new_msa = block(first_res, first_score)
        else:
            break
    print(first_score)
    for r in first_res:
        print(r)

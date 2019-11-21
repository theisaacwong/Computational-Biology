import sys
import copy
import math
import argparse
import time
import os


parser = argparse.ArgumentParser()
parser.add_argument("distance_matrix_path", help="path to the distance matrix file")
parser.add_argument("cluster_coords_path", help="path to the cluster coordinates file")
parser.add_argument("-o", "--out", help="output file path", type=str, default="quantileMapping_out_" + str(int(time.time())) + ".txt")
parser.add_argument("--min", help="minimum length of arrays for filtering", type=int, default=3)
parser.add_argument("--max", help="maximum length of arrays for filtering", type=int, default=sys.maxsize)
parser.add_argument("--q1", help="first quantile in percent, or first index as integer", type=float, default=0.1)
parser.add_argument("--q2", help="second quantile in percent, or second index from end as integer", type=float, default=0.9)
parser.add_argument("--trim1", help="how many starting monomers to initially ignore in the array", type=int, default=0)
parser.add_argument("--trim2", help="how many ending monomers to initially ignore in the array", type=int, default=0)
args = parser.parse_args()

#  example args
#args = parser.parse_args(["dmau/1.688", "dmau_f_x_1.688.tsv", "--min", "4", "--q1", "1", "--q2", "1", "-o", "dmau_x_1.688_q1_q1_f4_03.txt"])
#args = parser.parse_args(["dmel_Xchrom_1.688_f.dist", "dmel_f_x_1.688.tsv", "--min", "3", "--q1", "2", "--q2", "2", "-o", "dmel_x_1.688_q2_q2_f3_01.txt"])
#args = parser.parse_args(["dsim_Xchrom_1.688_f.dist", "dsim_f_x_1.688.tsv", "--min", "3", "--q1", "2", "--q2", "2", "-o", "dsim_x_1.688_q2_q2_f3_01.txt"])
#args = parser.parse_args(["dsech_Xchrom_1.688_f.dist", "dsech_f_x_1.688.tsv", "--min", "3", "--q1", "2", "--q2", "2", "-o", "dsech_x_1.688_q2_q2_f3_01.txt"])


leadingTrim = args.trim1
laggingTrim = args.trim2
minFilterLength = args.min
maxFilterLength = args.max
firstQuantile = args.q1
lastQuantile = args.q2
distPath = args.distance_matrix_path
clustPath = args.cluster_coords_path
outPath = args.out

distMap = {}
clustMap = {}

index1 = 0
index2 = 0
index3 = 0
index4 = 0

# old
# with open(distPath) as distFile:  # open file
#     headerLinee = distFile.readline().strip().split("\t")  # grab header line, by strip()-ing, first tab character is removed
#     for line in distFile:  # for the remaining lines
#         linee = line.strip().split("\t")  # split the line into an array by tabs
#         for i in range(0, len(headerLinee)):  # for each coord in the header array
#             start = int(linee[0].split("-")[1])
#             end = int(headerLinee[i].split("-")[0].split(":")[1])
#             physDist = str(abs(end - start))
#             genDist = str(linee[i+1])
#             distMap[headerLinee[i] + " vs " + linee[0]] = genDist + "\t" + physDist  # key: header + linee[0]; value: distance
temp1 = 0
ultimateHeaderLinee = []
# new, need to make sure parsing is done right
for distMatrix in os.listdir(distPath):
    with open(distPath + "/" + distMatrix) as distFile:  # open file
        headerLinee = distFile.readline().strip().split("\t")  # grab header line, by strip()-ing, first tab character is removed
        ultimateHeaderLinee = ultimateHeaderLinee + headerLinee
        for line in distFile:  # for the remaining lines
            linee = line.strip().split("\t")  # split the line into an array by tabs
            secondCoordTemp1 = linee[0].split("X")[1]  # parse the coordinates
            secondCoordStart = secondCoordTemp1.split("-")[0].split(":")[1]
            secondCoordEnd = secondCoordTemp1.split("-")[1]
            secondCoord = "X:" + secondCoordStart + "-" + secondCoordEnd
            for i in range(0, len(headerLinee)):  # for each coord in the header array
                headerCoordTemp1 = headerLinee[i].split("X")[1]
                headerCoordStart = headerCoordTemp1.split("-")[0].split(":")[1]
                headerCoordEnd = headerCoordTemp1.split("-")[1]
                headerCoord = "X:" + headerCoordStart + "-" + headerCoordEnd

                start = int(headerCoordStart)
                end = int(headerCoordEnd)
                physDist = str(abs(end - start))
                genDist = str(linee[i+1])
                if(genDist == "NA"):
                    genDist = "nan"
                distMap[headerCoord + " vs " + secondCoord] = genDist + "\t" + physDist  # key: header + linee[0]; value: distance


                if headerCoord == secondCoord:
                    temp1 += 1

headerLinee = ultimateHeaderLinee
# end new

for k, v in distMap.items():
    print(k, v)
print(temp1)
print(len(distMap.keys()))

with open(clustPath) as clustFile:
    for i, line in enumerate(clustFile):
        clustMap[i] = line.strip().split("\t")

# filter the clusters, done separately from parsing on purpose
for i in range(0, len(clustMap.keys())):
    if int(clustMap[i][1]) < minFilterLength or int(clustMap[i][1]) > maxFilterLength:
        del clustMap[i]

# Pseudo code for the following portion
# for each cluster
#   get all individual monomers
#   all-by-all comparison, generate a square graph
#   separate square graph into quantiles
#   for each quantile, average genetic distance and physical distance
#   return one square graph
# for all individual square graphs
#   average genetic distance
#   weight by physical distance?????

listOf3x3 = []
sizeArray = []
bandArray = []

# for storing the coordinates of each array
monomer_first_start = []
monomer_first_end = []
monomer_last_start = []
monomer_last_end = []
monomer_second_first = []
monomer_second_last = []


for k in clustMap.keys():
    # first need to generate the smaller 'window' square
    value3x3 = [[0 for x in range(3)] for y in range(3)]
    # get the size of the current cluster
    size = int(clustMap[k][1])
    # make 2D arrays to store the key/value strings, key: X:__ vs X:__   value: genetic distance
    key2DArray = [[" " for x in range(size)] for y in range(size)]
    value2DArray = [[" " for x in range(size)] for y in range(size)]
    # currentCoords, all the coordinates in between the first and last monomers
    currCoords = []
    # coordinates of the first monomer
    first_start = clustMap[k][7]
    first_end = clustMap[k][8]

    # get the band that the array is from
    bandArray.append(str(clustMap[k][11]))

    # re-create the header coordinates of the first monomer
    # old
    # startIndex = headerLinee.index("X:" + first_start + "-" + first_end)
    startIndex = -1
    for i in range(len(headerLinee)):
        if ("X:" + first_start + "-" + first_end) in headerLinee[i]:
            startIndex = i

    if startIndex == -1:
        print("error")
        print("X:" + first_start + "-" + first_end)

    # get all the coordinates in between the first and last monomers by subsetting from the header array
    for i in range(size):
        currCoords.append("X" + headerLinee[startIndex + i].split("X")[1])

    # trim the currCoords array to only look at monomers of interest
    tempCoords = []
    size = size - leadingTrim - laggingTrim
    for i in range(size):
        tempCoords.append(currCoords[i + leadingTrim])
    currCoords = tempCoords
    monomer_first_start.append(tempCoords[0])
    monomer_last_start.append(tempCoords[size-1])
    monomer_second_first.append(tempCoords[1])
    monomer_second_last.append(tempCoords[size-2])

    # re-create all the keys and get their values for each pair of coordinates
    for i in range(size):
        for j in range(size):
            key2DArray[i][j] = currCoords[i] + " vs " + currCoords[j]
            value2DArray[i][j] = float(distMap[key2DArray[i][j]].split("\t")[0])

    # now we need to separate into quantiles, and average

    value3x3 = [[0 for x in range(3)] for y in range(3)]

    # check to see if the user wants quantiles defined by percent or by number
    if firstQuantile >= 1:
        firstQuantileIndex = int(firstQuantile)
    else:
        firstQuantileIndex = max(int(math.ceil(float(size) * firstQuantile)), 1)  # rounding up

    if lastQuantile >= 1:
        lastQuantileIndex = size - int(lastQuantile)
    else:
        lastQuantileIndex = min(int(math.ceil(float(size) * lastQuantile)), size-1)  # rounding down

    # I am really really really really sorry for how dumb this is
    # I know I could slice the original 2D matrix with sum([a:b, x:y])
    # but you overestimate my intelligence
    # might've been easier in R
    total, counter = 0.0, 0.0  # first row first column - top left
    print(size, firstQuantileIndex, lastQuantileIndex)
    for i in range(leadingTrim, firstQuantileIndex):
        for j in range(0, firstQuantileIndex):
            if i != j:  # don't include the diagonal
                total += value2DArray[i][j]
                counter += 1
    if counter == 0 and total == 0:
        value3x3[0][0] = 0
    else:
        value3x3[0][0] = total / counter

    total, counter = 0.0, 0.0  # first row second column - top middle
    for i in range(0, firstQuantileIndex):
        for j in range(firstQuantileIndex, lastQuantileIndex):
            if i != j:  # don't include the diagonal
                total += value2DArray[i][j]
                counter += 1
    if counter == 0 and total == 0:
        value3x3[0][1] = 0
    else:
        value3x3[0][1] = total / counter

    total, counter = 0.0, 0.0  # first row third column - top right
    for i in range(0, firstQuantileIndex):
        for j in range(lastQuantileIndex, size):
            if i != j:  # don't include the diagonal
                total += value2DArray[i][j]
                counter += 1
    if counter == 0 and total == 0:
        value3x3[0][2] = 0
    else:
        value3x3[0][2] = total / counter

    total, counter = 0.0, 0.0  # second row first column - middle left
    for i in range(firstQuantileIndex, lastQuantileIndex):
        for j in range(0, firstQuantileIndex):
            if i != j:  # don't include the diagonal
                total += value2DArray[i][j]
                counter += 1
    if counter == 0 and total == 0:
        value3x3[1][0] = 0
    else:
        value3x3[1][0] = total / counter

    total, counter = 0.0, 0.0  # second row second column - middle middle
    for i in range(firstQuantileIndex, lastQuantileIndex):
        for j in range(firstQuantileIndex, lastQuantileIndex):
            if i != j:  # don't include the diagonal
                total += value2DArray[i][j]
                counter += 1
    if counter == 0 and total == 0:
        value3x3[1][1] = 0
    else:
        value3x3[1][1] = total / counter

    total, counter = 0.0, 0.0  # second row third column - middle right
    for i in range(firstQuantileIndex, lastQuantileIndex):
        for j in range(lastQuantileIndex, size):
            if i != j:  # don't include the diagonal
                total += value2DArray[i][j]
                counter += 1
    if counter == 0 and total == 0:
        value3x3[1][2] = 0
    else:
        value3x3[1][2] = total / counter

    total, counter = 0.0, 0.0  # third row first column - bottom left
    for i in range(lastQuantileIndex, size):
        for j in range(0, firstQuantileIndex):
            if i != j:  # don't include the diagonal
                total += value2DArray[i][j]
                counter += 1
    if counter == 0 and total == 0:
        value3x3[2][0] = 0
    else:
        value3x3[2][0] = total / counter

    total, counter = 0.0, 0.0  # third row second column - bottom middle
    for i in range(lastQuantileIndex, size):
        for j in range(firstQuantileIndex, lastQuantileIndex):
            if i != j:  # don't include the diagonal
                total += value2DArray[i][j]
                counter += 1
    if counter == 0 and total == 0:
        value3x3[2][1] = 0
    else:
        value3x3[2][1] = total / counter

    total, counter = 0.0, 0.0  # third row third column - bottom right
    for i in range(lastQuantileIndex, size):
        for j in range(lastQuantileIndex, size):
            if i != j:  # don't include the diagonal
                total += value2DArray[i][j]
                counter += 1
    if counter == 0 and total == 0:
        value3x3[2][2] = 0
    else:
        value3x3[2][2] = total / counter

    listOf3x3.append(copy.deepcopy(value3x3))  # save this 3x3 array - add it to the list of 3x3 arrays
    sizeArray.append(size)

meanSize = sum(sizeArray) / float(len(sizeArray))
# get the mean of all the 3x3 plots
mean3x3 = [[0 for x in range(3)] for y in range(3)]
numArrays = len(listOf3x3)
for i in range(0, numArrays):
    for j in range(0, len(listOf3x3[i])):
        for k in range(0, len(listOf3x3[i][j])):
            mean3x3[j][k] += listOf3x3[i][j][k]

for i in range(0, len(mean3x3)):
    for j in range(0, len(mean3x3[i])):
        mean3x3[i][j] /= numArrays

for row in mean3x3:
    # Loop over columns.
    for column in row:
        print(column, end=" \t")
    print(end="\n")

headers = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
headers = ["first-first", "first-middle", "first-last",
           "middle-first", "middle-middle", "middle-last",
           "last-first", "last-middle", "last-last", "size", "band",
           "first-start", "second-start", "last-start", "last-second"]

# might've been easier to do in R
with open(outPath, 'w') as out:
    out.write("\t".join(headers))
    out.write("\n")
    for row in mean3x3:
        out.write("\t".join([str(x) for x in row]) + "\t")
    out.write(str(meanSize) + "\t" + "NA" + "\tNA\tNA\tNA\tNA")
    out.write("\n")
    for i in range(0, numArrays):
        for j in range(0, len(listOf3x3[i])):
            for k in range(0, len(listOf3x3[i][j])):
                out.write(str(listOf3x3[i][j][k]))
                if j != len(listOf3x3[i]) and k != len(listOf3x3[i][j]):
                    out.write("\t")
        out.write(str(sizeArray[i]) + "\t" + bandArray[i] + "\t")
        out.write(str(monomer_first_start[i]) + "\t")
        out.write(str(monomer_last_start[i]) + "\t")
        out.write(str(monomer_second_first[i]) + "\t")
        out.write(str(monomer_second_last[i]))
        out.write("\n")


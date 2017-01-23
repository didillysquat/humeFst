from django.shortcuts import render
from django.http import HttpResponse
from .models import its2Seq, dataSubmission
from django.db.models import Q
import json # This is what you use to json a dictionary object
from django.core import serializers # this is what you use to json a query set



def home(request):
    methodtype = request.method
    return render(request, 'dbApp/home.html', {'method':methodtype})

def learning(request):

    return render(request, 'dbApp/learning.html')

def summary(request):
    methodtype = request.method
    return render(request, 'dbApp/summary.html', {'method':methodtype})


def bioDivAnalysis(request):
    return render(request, 'dbApp/bioDivAnalysis.html', {'checker': 'This is the bioDivAnalysis page'})

def dataInputValidation(request):
    a = None
    # To start with let's see if we can access a text document and maybe print out a line from it to a file.
    if request.method == 'POST':
        if request.POST['dataType'] == 'fastQ':
            # Then here we read in one file
            instance = dataSubmission(submissionType='fastQ', submittingUser='Bob', submittedDataRaw=request.FILES['fastQFile'])
            instance.save()
            currentDataSubmissionID = instance.id
            letsSeeFile = instance.submittedDataRaw
            a = 'something'
        elif request.POST['dataTYpe'] == 'fasta':
            # Then here we read in multiple files
            c = 'something else'
        a = 5
    return render(request, 'dbApp/biodDivAnalysis.html', {'something': a, 'somethingElse': a})

def seqQuery(request):

    # Initialize querySetToReturn as None
    querySetToReturn = None
    queryDictPrint = {}


    if request.method == 'GET':
        # Then this is the first vist to the page
        # Return all sequences into the table
        querySetToReturn = its2Seq.objects.all().order_by('-occurrence') # This returns a queryset that can then be refined/filtered
        queryDictPrint = {'seq_namePrint': 'All', 'seq_cladePrint': 'All',
                          'seq_lengthPrint': 'All',
                          'seq_occurrencePrint': 'All'}


    elif request.method == 'POST':
        # Then we get the query parameters and we give back a seqCollection according to these parameters
        # First check to see which fields have a value and initialize query arguments accordingly
        queryDictPrint = {'seq_namePrint': request.POST['seq_name'], 'seq_cladePrint': request.POST['seq_clade'], 'seq_lengthPrint': request.POST['seq_length'], 'seq_occurrencePrint': request.POST['seq_occurrence']}
        queryDict = {}
        seq_name_reply = request.POST['seq_name']
        if seq_name_reply == '':
            queryDict['seq_name'] = 'ALL'
        else:
            seqNameEx = seq_name_reply
            # Here we need to extract each of the seqnames if it is a list e.g. C3, D1, A1
            queryDict['seq_name'] = [x.strip() for x in seqNameEx.split(',')]


        seq_clade_reply = request.POST['seq_clade']
        if seq_clade_reply == 'ALL':
            queryDict['seq_clade'] = 'ALL'
        else:
            queryDict['seq_clade'] = seq_clade_reply


        seq_length_reply = request.POST['seq_length']
        if seq_length_reply == '':
            queryDict['seq_length'] = 'ALL'
        else:
            if '-' in seq_length_reply:
                queryDict['seq_length'] = seq_length_reply.split('-')
            else:
                queryDict['seq_length'] = seq_length_reply


        seq_occurrence_reply = request.POST['seq_occurrence']
        if seq_occurrence_reply == '':
            queryDict['seq_occurrence'] = 'ALL'
        else:
            if '-' in seq_occurrence_reply:
                queryDict['seq_occurrence'] = seq_occurrence_reply.split('-')
            else:
                queryDict['seq_occurrence'] = seq_occurrence_reply

        seq_sortby_reply = request.POST['seq_sortby']
        if seq_sortby_reply == 'Sequence occurrence':
            queryDict['seq_sortby'] = '-occurrence'
        elif seq_sortby_reply == 'Sequence name':
            queryDict['seq_sortby'] = 'name'
        elif seq_sortby_reply == 'Sequence length':
            queryDict['seq_sortby'] = '-length'
        elif seq_sortby_reply == 'Sequence clade':
            queryDict['seq_sortby'] = 'clade'


        # Here we now need to make the queryset as defined by the terms in queryDict
        querySetToReturn = its2Seq.objects.all()

        # Filter by name
        if queryDict['seq_name'] != 'ALL':
            predicates = [('name', seqName) for seqName in queryDict['seq_name']]
            qList = [Q(x) for x in predicates]
            args = Q()
            for each_args in qList:
                args = args | each_args
            querySetToReturn = querySetToReturn.filter(*(args,))

        # Filter by clade
        if queryDict['seq_clade'] != 'ALL':
            querySetToReturn = querySetToReturn.filter(clade=queryDict['seq_clade'])

        # Filter by length
        if queryDict['seq_length'] != 'ALL':
            if len(queryDict['seq_length']) > 1:
                minValue = queryDict['seq_length'][0]
                maxValue = queryDict['seq_length'][1]
                # Then this is a range, min first then max
                querySetToReturn = querySetToReturn.filter(length__gte=minValue, length__lte=maxValue)
            else:
                # Then this is a single value and we need to work out the operator
                queryOperator = queryDict['seq_length'][:2]
                if queryOperator == '>=':
                    querySetToReturn = querySetToReturn.filter(length__gte=queryDict['seq_length'][2:])
                elif queryOperator == '<=':
                    querySetToReturn = querySetToReturn.filter(length__lte=queryDict['seq_length'][2:])
                elif queryOperator[0] == '>':
                    querySetToReturn = querySetToReturn.filter(length__gt=queryDict['seq_length'][1:])
                elif queryOperator[0] == '<':
                    querySetToReturn = querySetToReturn.filter(length__lt=queryDict['seq_length'][1:])

        # Filter by occurrence
        if queryDict['seq_occurrence'] != 'ALL':
            if len(queryDict['seq_occurrence']) > 1:
                minValue = queryDict['seq_occurrence'][0]
                maxValue = queryDict['seq_occurrence'][1]
                # Then this is a range, min first then max
                querySetToReturn = querySetToReturn.filter(occurrence__gte=minValue, occurrence__lte=maxValue)
            else:
                # Then this is a single value and we need to work out the operator
                queryOperator = queryDict['seq_occurrence'][:2]
                if queryOperator == '>=':
                    querySetToReturn = querySetToReturn.filter(occurrence__gte=queryDict['seq_occurrence'][2:])
                elif queryOperator == '<=':
                    querySetToReturn = querySetToReturn.filter(occurrence__lte=queryDict['seq_occurrence'][2:])
                elif queryOperator[0] == '>':
                    querySetToReturn = querySetToReturn.filter(occurrence__gt=queryDict['seq_occurrence'][1:])
                elif queryOperator[0] == '<':
                    querySetToReturn = querySetToReturn.filter(occurrence__lt=queryDict['seq_occurrence'][1:])

        # Sort according to seq_sortby

        querySetToReturn = querySetToReturn.order_by('{0}'.format(queryDict['seq_sortby']))


    return render(request, 'dbApp/seqQuery.html', {'querySet':querySetToReturn, 'queryDictPrint':queryDictPrint})



def seqModal(request):
    seq_id = None
    context = 'Apple'
    if request.method == 'GET':
        sequenceID = request.GET['seq_id']

    if sequenceID:
        seq = its2Seq.objects.get(id=int(sequenceID))
        title = 'Sequence: {0}'.format(seq.name)
        dataReturn = {'title': title, 'sequence': seq.sequence, 'fastaTitle': '>' + seq.name}
        dataReturn = json.dumps(dataReturn)
    return HttpResponse(dataReturn)



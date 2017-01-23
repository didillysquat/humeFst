from django.db import models

class its2Seq(models.Model):
    name = models.CharField(max_length=20)
    sequence = models.TextField()
    length = models.IntegerField()
    clade = models.CharField(max_length=10)
    occurrence = models.IntegerField()

    def __str__(self):
        return self.name

class dataSubmission(models.Model):

    submissionType_choices = (('fasta', 'fasta'),('fastQ', 'fastQ'))
    submissionType = models.CharField(max_length=5, choices=submissionType_choices)
    submittingUser = models.CharField(max_length=30, default='Bob')
    submittedDataRaw = models.FileField()


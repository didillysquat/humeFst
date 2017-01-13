from django.db import models

class its2Seq(models.Model):
    name = models.CharField(max_length=20)
    sequence = models.TextField()
    length = models.IntegerField()
    clade = models.CharField(max_length=10)
    occurrence = models.IntegerField()

    def __str__(self):
        return self.name

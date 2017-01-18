from django.conf.urls import url
from django.views.generic import ListView, DetailView
from .models import its2Seq
from . import views

app_name = 'dbApp'

urlpatterns = [
    url(r'^$', views.home, name='home'),
    url(r'^seqQuery$', views.seqQuery, name='seqQuery'),
    url(r'^summary', views.summary, name='summary'),
    url(r'^(?P<pk>\d+)$', DetailView.as_view(model=its2Seq, template_name = 'dbApp/its2SeqDetail.html')),
    url(r'^seqModal/$', views.seqModal, name='seqModal'),
    url(r'^bioDivAnalysis$', views.bioDivAnalysis, name='bioDivAnalysis')
]




# url(r'^$', ListView.as_view(queryset=its2Seq.objects.all().order_by("-occurance")[:25], template_name="dbApp/seqList.html")),
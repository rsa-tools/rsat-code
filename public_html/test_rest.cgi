use Mojolicious::Lite;

get '/questions/(:question_id)' => sub {
    my $self = shift;
    my $result = {};
    # do stuff with $result based on $self->stash('question_id')
    result->{answer} = "Hello Mojolicious";
    return $self->render_json($result)
}

app->start;

